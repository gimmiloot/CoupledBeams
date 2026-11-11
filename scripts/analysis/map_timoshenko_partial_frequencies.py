"""Partial Timoshenko CS/CC frequencies with immutable coupled-map inputs.

New workflow contract: recover certified CC catalogs, solve only missing CS
lengths, then plot saved tables. Reuses the existing finite-rod library through
a boundary adapter. No coupled spectrum, FEM, Ritz or tracking is called.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
DEFAULT=ROOT/'data/input/timoshenko_partial_m15_0.json'


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path): return json.loads(Path(path).read_text(encoding='utf-8'))
def write(path,obj):
    def scalar(value):
        if isinstance(value,np.generic): return value.item()
        if isinstance(value,np.ndarray): return value.tolist()
        raise TypeError(type(value).__name__)
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False,default=scalar)+'\n',encoding='utf-8')
def rows(path):
    with Path(path).open(encoding='utf-8',newline='') as f: return list(csv.DictReader(f))
def table(path,items):
    with Path(path).open('w',encoding='utf-8',newline='') as f:
        if items:
            w=csv.DictWriter(f,fieldnames=list(items[0]));w.writeheader();w.writerows(items)


def compute(config_path,output):
    from scripts.lib import timoshenko_partial_spectra as part
    from scripts.lib.rl_full_rr_aspect import make_model
    from scripts.lib import mindlin_herrmann_longitudinal as mh
    config=read(config_path);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    assert config['material']==dict(epsilon=.005,aspect=1.,E=1.,rho=1.,l=1.,nu=.3,ks=5/6,variant='RL_FULL')
    assert config['source_series']=='M15_0' and config['Lambda_top']==11.
    (output/'points').mkdir(exist_ok=True)
    previous=read(output/'manifest.json') if (output/'manifest.json').exists() else {}
    previous_records={r['record']:r for r in previous.get('points',[])} if previous.get('status')=='PASS' else {}
    policy=read(ROOT/config['policy_source'])['policy']
    model=make_model(**config['material']).arms[0].carrier
    p=model.coefficients;section=model.section
    omega_top=config['material']['epsilon']*config['Lambda_top']**2
    inventory=read(ROOT/config['catalog_inventory']); expected=read(ROOT/config['source_hashes'])
    lengths=sorted(round(i['length'],12) for i in inventory)
    assert lengths==[round(.1+.005*i,12) for i in range(361)]
    contract=dict(version=part.VERSION,material=config['material'],policy=policy,Lambda_top=config['Lambda_top'],omega_top=omega_top,
        state_order=['W','Psi','Q','M'],section=dict(b=section.width,h=section.thickness,A=section.area,I=section.inertia),
        coefficients={k:p[k] for k in ('B','S','m','r')},normalization='omega=epsilon*Lambda**2',
        adapter_sha256=sha(part.__file__),basis_sha256=sha(mh.__file__),inventory_sha256=sha(ROOT/config['catalog_inventory']))
    identity=hashlib.sha256(json.dumps(contract,sort_keys=True).encode()).hexdigest()
    source_hashes={str(Path(config_path).relative_to(ROOT)):sha(config_path),config['catalog_inventory']:sha(ROOT/config['catalog_inventory']),
                   config['source_hashes']:sha(ROOT/config['source_hashes']),config['policy_source']:sha(ROOT/config['policy_source'])}
    records=[]; frequency_rows=[];checks=[];counts=[];new_cs=reused=0
    for index,item in enumerate(sorted(inventory,key=lambda i:i['length'])):
        L=round(item['length'],12); source=ROOT/item['record']
        assert sha(source)==expected[item['record']]
        source_hashes[item['record']]=sha(source)
        original=read(source)
        cc=next(c for c in original['catalogs'].values() if c['catalog_id']==item['catalog_id'])
        assert cc['status']=='PASS' and cc['upper']>=omega_top and not cc['search']['failed_intervals']
        assert len(cc['poles'])==cc['count_bound']['upper_count']==len(cc['roots'])
        np.testing.assert_allclose(cc['arm_contract']['length'],L,rtol=0,atol=2e-15)
        for a,b in [('EI','B'),('KGA','S'),('rhoA','m'),('rhoI','r')]:
            np.testing.assert_allclose(cc['arm_contract'][a],p[b],rtol=2e-14,atol=0)
        for key,val in config['material'].items(): assert cc['arm_contract']['material'][key]==val
        assert mh.finite_count_upper_bound(model,L,.001*np.pi,'timoshenko')['upper_count']==0
        for bc in ('CC','CS'):
            name=f'{bc}_L{L:.3f}'; path=output/'points'/(name+'.json')
            if path.exists():
                if previous_records:
                    assert sha(path)==previous_records[path.relative_to(output).as_posix()]['sha256']
                rec=read(path)
                assert rec['status']=='PASS' and rec['identity']==identity and rec['source_sha256']==sha(source)
                assert rec['bc']==bc and rec['length']==L and rec['catalog_id']==item['catalog_id']
                assert sha(output/rec['npz'])==rec['npz_sha256']
                reused+=1
            else:
                try:
                    if bc=='CC':
                        roots=[q for q in cc['roots'] if q['omega']<=omega_top]
                        certificate=dict(status='PASS',method='immutable certified CC catalogue',catalog_id=item['catalog_id'],
                            upper=cc['upper'],upper_count=cc['count_bound']['upper_count'],search=cc['search'],new_root_searches=0)
                    else:
                        roots,certificate=part.solve_cs(model,L,cc,omega_top,policy)
                        new_cs+=certificate['search_calls']
                    modes,quality=part.verify_modes(model,L,roots,bc,policy)
                    control=None
                    if bc=='CS' and L in config['control_lengths'] and roots:
                        finer=dict(policy,pole_scan_intervals=2*policy['pole_scan_intervals'])
                        refined,_=part.solve_cs(model,L,cc,omega_top,finer)
                        error=max(abs(a['omega']/b['omega']-1) for a,b in zip(roots,refined))
                        assert len(roots)==len(refined) and error<=policy['same_physics_frequency_tol']
                        control=dict(scan_intervals=[policy['pole_scan_intervals'],finer['pole_scan_intervals']],relative_frequency_change=error)
                    arrays={f'mode{k}_{key}':m[key] for k,m in enumerate(modes,1) for key in ('points','weights','values','gradients','ends','coefficients')}
                    npz=output/'points'/(name+'.npz');np.savez_compressed(npz,**arrays)
                    rec=dict(status='PASS',identity=identity,length=L,bc=bc,roots=roots,certificate=certificate,quality=quality,
                        control=control,source=item['record'],source_sha256=sha(source),catalog_id=item['catalog_id'],
                        npz=npz.relative_to(output).as_posix(),npz_sha256=sha(npz))
                    write(path,rec)
                except Exception as exc:
                    fail=dict(status='FAIL',length=L,bc=bc,identity=identity,reason=str(exc),completed=len(records))
                    write(output/(name+'_FAIL.json'),fail);write(output/'manifest.json',fail);raise
            records.append(dict(length=L,bc=bc,record=path.relative_to(output).as_posix(),sha256=sha(path),npz=rec['npz'],npz_sha256=rec['npz_sha256']))
            inband=[r for r in rec['roots'] if r['omega']<=omega_top]
            counts.append(dict(length=L,bc=bc,status='PASS',roots_in_band=len(inband),roots_saved=len(rec['roots']),
                gram_max_error=rec['quality']['gram_max_error'],source=item['record'],catalog_id=item['catalog_id']))
            for k,r in enumerate(rec['roots'],1):
                frequency_rows.append(dict(length=L,bc=bc,rank=k,omega=r['omega'],Lambda=float(np.sqrt(r['omega']/.005)),
                    role='plot' if r['omega']<=omega_top else 'guard',record=path.relative_to(output).as_posix(),
                    source=item['record'],source_sha256=rec['source_sha256'],catalog_id=item['catalog_id']))
            for d in rec['quality']['modes']:
                checks.append(dict(length=L,bc=bc,**{k:v for k,v in d.items() if k!='gates'},gates_pass=all(d['gates'].values())))
        if index%30==0: print(f'{index+1}/361 lengths PASS; L={L:.3f}',flush=True)
    # Monotonicity of each isolated-beam sorted eigenfrequency, not tracking.
    for bc in ('CS','CC'):
        for rank in range(1,8):
            seq=sorted((r for r in frequency_rows if r['bc']==bc and r['rank']==rank),key=lambda r:r['length'])
            assert all(b['omega']<a['omega'] for a,b in zip(seq,seq[1:]))
    table(output/'partial_frequencies.csv',frequency_rows);table(output/'length_status.csv',counts);table(output/'quality.csv',checks)
    cmap=ROOT/config['source_map']; cm=read(cmap/'manifest.json');source_hashes[(cmap/'manifest.json').relative_to(ROOT).as_posix()]=sha(cmap/'manifest.json')
    assert cm['status']=='PASS' and sha(cmap/'spectrum.csv')==cm['artifacts']['spectrum.csv']
    coupled=[r for r in rows(cmap/'spectrum.csv') if r['map']==config['source_series']]
    assert len(coupled)==181*7
    for point in (r for r in cm['points'] if r['map']==config['source_series']):
        rp=ROOT/point['record'];assert sha(rp)==point['record_sha256'];r=read(rp)
        assert r['status']=='PASS' and len(r['roots'])==7 and not r['search']['unresolved_intervals']
        assert (r['beta_deg'],r['tau'],r['contract']['K_plot'],r['contract']['K_guard'])==(15.,0.,6,7)
        assert all(r['contract']['model'][k]==v for k,v in config['material'].items())
        for row in (q for q in coupled if float(q['mu'])==r['mu']): assert float(row['Lambda'])==r['roots'][int(row['rank'])-1]['Lambda']
        npz=rp.parent.parent/r['shape_file'];assert sha(npz)==r['shape_sha256']
        source_hashes[point['record']]=sha(rp);source_hashes[npz.relative_to(ROOT).as_posix()]=sha(npz)
    table(output/'coupled_frequencies.csv',[r for r in coupled if int(r['rank'])<=6])
    lookup={(r['bc'],round(r['length'],12),r['rank']):r for r in frequency_rows}
    armrows=[]
    for i in range(181):
        mu=round(i*.005,12)
        for arm,L in [(1,round(1-mu,12)),(2,round(1+mu,12))]:
            for bc in ('CS','CC'):
                for rank in range(1,8):
                    r=lookup.get((bc,L,rank))
                    if r and r['role']=='plot': armrows.append(dict(mu=mu,arm=arm,**r))
    table(output/'arm_frequencies.csv',armrows)
    write(output/'contract.json',dict(identity=identity,physical=contract,config=config))
    write(output/'source_hashes.json',source_hashes)
    write(output/'code_provenance.json',{str(Path(f).relative_to(ROOT)):sha(f) for f in
        (__file__,part.__file__,mh.__file__,ROOT/'scripts/lib/rl_full_rr_aspect.py',ROOT/'scripts/lib/frequency_geometry_style.py')})
    manifest=dict(status='PASS',identity=identity,points=records,new_cs_searches_this_run=new_cs,reused_records_this_run=reused,
        cc_new_root_searches=0,coupled_new_solves=0,lengths=361,
        accepted_cs_search_geometries=sum(read(output/r['record'])['certificate']['search_calls'] for r in records if r['bc']=='CS'),
        roots_in_band={bc:sum(r['bc']==bc and r['role']=='plot' for r in frequency_rows) for bc in ('CS','CC')},
        cs_guards=sum(r['bc']=='CS' and r['role']=='guard' for r in frequency_rows),
        artifacts={p.relative_to(output).as_posix():sha(p) for p in output.rglob('*') if p.is_file() and p.name!='manifest.json' and 'figures' not in p.parts})
    manifest['resolved_technical_events']=[dict(file=p.name,reason=read(p)['reason'],
        resolution='NumPy JSON scalar serialization fixed; completed records reused; physical gates unchanged')
        for p in output.glob('*_FAIL.json') if 'JSON serializable' in read(p)['reason']]
    write(output/'manifest.json',manifest)
    print(json.dumps({k:v for k,v in manifest.items() if k not in ('points','artifacts')}),flush=True)
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=DEFAULT)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--plot-only',action='store_true')
    parser.add_argument('--figure-output',type=Path)
    args=parser.parse_args()
    if args.plot_only:
        render(args.output,args.figure_output or args.output/'figures')
    else:
        compute(args.config,args.output)


def render(output,figures):
    # This route reads completed tables. It does not import the partial solver.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from scripts.lib import frequency_geometry_style as style
    output,figures=Path(output),Path(figures)
    manifest=read(output/'manifest.json');assert manifest['status']=='PASS'
    for name in ('contract.json','partial_frequencies.csv','coupled_frequencies.csv','arm_frequencies.csv'):
        assert sha(output/name)==manifest['artifacts'][name]
    if figures.exists() and any(figures.iterdir()): raise FileExistsError('Choose an empty new figure directory')
    figures.mkdir(parents=True,exist_ok=True)
    compare_saved(output,figures)
    coupled=rows(output/'coupled_frequencies.csv');arms=rows(output/'arm_frequencies.csv')
    contract=read(output/'contract.json');mus=sorted({float(r['mu']) for r in coupled})
    material=contract['physical']['material']
    geometry=style.FrequencyGeometry(mus,axis='MU',fixed=dict(beta_deg=15.,tau=0.),material=material)
    plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','DejaVu Serif'],
        'mathtext.fontset':'stix','font.size':11,'axes.linewidth':.65,'pdf.fonttype':42})
    colors={1:'#766052',2:'#71817b'}
    for variant,families in [('A',['CS']),('B',['CS','CC'])]:
        fig,ax=plt.subplots(figsize=(8.1,5.1))
        for bc in families:
            for arm in (1,2):
                for rank in range(1,8):
                    points={float(r['mu']):float(r['Lambda']) for r in arms if r['bc']==bc and int(r['arm'])==arm and int(r['rank'])==rank}
                    if not points: continue
                    values=[points.get(mu,np.nan) for mu in mus]
                    style.plot_partial_frequency(ax,mus,values,color=colors[arm],label='_nolegend_',lw=.85,
                        marker='s' if bc=='CC' else None,markersize=2.4,markerfacecolor='white',markeredgewidth=.55,markevery=11,zorder=1)
        for rank in range(1,7):
            q=sorted((r for r in coupled if int(r['rank'])==rank),key=lambda r:float(r['mu']))
            values=[float(r['Lambda']) for r in q]
            style.plot_frequency_curve(ax,geometry,values,color='#173e5c',lw=1.35,zorder=3)
            ax.text(.907,values[-1],str(rank),fontsize=10,va='center',color='#173e5c',clip_on=False)
        handles=[Line2D([],[],color='#173e5c',lw=1.35,label='Конструкция: шесть частот')]
        for bc in families:
            for arm in (1,2):
                handles.append(Line2D([],[],color=colors[arm],ls=':',lw=.9,
                    marker='s' if bc=='CC' else None,markersize=3.5,markerfacecolor='white',
                    label=f'{bc}, стержень {arm}'))
        ax.legend(handles=handles,loc='lower left',bbox_to_anchor=(0.,1.01),ncol=3,frameon=False,fontsize=9,
                  handlelength=2.8,columnspacing=1.4,borderaxespad=0.)
        ax.set(xlim=(0,.9),ylim=(2,11),xlabel=r'$\mu$',ylabel=r'$\Lambda$')
        ax.set_xticks(np.arange(0,1,.1));ax.set_yticks(range(2,12));ax.tick_params(direction='out',width=.65,length=3)
        ax.spines[['top','right']].set_visible(False)
        fig.subplots_adjust(left=.075,right=.955,bottom=.105,top=.85 if variant=='B' else .90)
        caption=style.add_geometry_caption(fig)
        name=f'frequency_{variant}_'+('CS' if variant=='A' else 'CS_CC')
        fig.savefig(figures/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None})
        fig.savefig(figures/(name+'.png'),dpi=220,bbox_inches='tight');plt.close(fig)
    style.write_geometry_diagnostics(figures,{'M15_0':geometry})
    write(figures/'manifest.json',dict(status='PASS',calculation_mode='plot_only',new_spectral_calls=0,
        spectrum_semantics='sorted_positions; common color; labels 1-6 are current order, not descendants',
        coupled_positions=[1,2,3,4,5,6],partial_families=['CS1','CS2','CC1','CC2'],
        data_sources={name:sha(output/name) for name in ('partial_frequencies.csv','arm_frequencies.csv','coupled_frequencies.csv')},
        partial_style='muted arm colors, dotted CS; dotted CC with open squares; never main dashed style',
        geometric_caption=caption,outputs={p.name:sha(p) for p in figures.iterdir() if p.is_file()}))
    print('Saved two PDF/PNG variants in '+str(figures))


def compare_saved(output,destination):
    """Crossing brackets from the existing length grid, not new eigenvalues."""
    arms=rows(Path(output)/'arm_frequencies.csv');coupled=rows(Path(output)/'coupled_frequencies.csv')
    intersections=[]
    for bc in ('CS','CC'):
        for n in range(1,7):
            for m in range(n+1,7):
                a={float(r['mu']):float(r['Lambda']) for r in arms if r['bc']==bc and int(r['arm'])==1 and int(r['rank'])==n}
                b={float(r['mu']):float(r['Lambda']) for r in arms if r['bc']==bc and int(r['arm'])==2 and int(r['rank'])==m}
                grid=sorted(set(a)&set(b))
                for left,right in zip(grid,grid[1:]):
                    dl,dr=a[left]-b[left],a[right]-b[right]
                    if dl*dr<0:
                        t=-dl/(dr-dl)
                        intersections.append(dict(bc=bc,rank_arm1=n,rank_arm2=m,mu_left=left,mu_right=right,
                            difference_left=dl,difference_right=dr,mu_linear_estimate=left+t*(right-left),
                            Lambda_linear_estimate=a[left]+t*(a[right]-a[left]),method='intersection of saved line segments; not an off-grid spectral solution'))
    table(Path(destination)/'partial_intersections.csv',intersections)
    mus=sorted({float(r['mu']) for r in coupled})
    f=np.array([[float(next(r['Lambda'] for r in coupled if float(r['mu'])==mu and int(r['rank'])==k)) for k in range(1,7)] for mu in mus])
    gaps=[];turns=[]
    for k in range(5):
        d=f[:,k+1]-f[:,k];i=int(np.argmin(d))
        gaps.append(dict(positions=[k+1,k+2],mu=mus[i],Lambda1=f[i,k],Lambda2=f[i,k+1],gap=d[i]))
    for k in range(6):
        for i in range(1,len(mus)-1):
            if (f[i,k]-f[i-1,k])*(f[i+1,k]-f[i,k])<0:
                turns.append(dict(position=k+1,mu=mus[i],Lambda=f[i,k],kind='max' if f[i,k]>f[i-1,k] else 'min'))
    write(Path(destination)/'comparison.json',dict(intersections=intersections,grid_gap_minima=gaps,grid_extrema=turns,
        interpolation_is_not_a_spectral_solution=True,coupled_solves=0,tracking=False))


if __name__=='__main__': main()
