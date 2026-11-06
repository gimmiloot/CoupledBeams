"""HS-B pure plotting from accepted CSV: no imports of spectral solvers."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def render(folder):
    folder=Path(folder)
    def table(name):
        with (folder/name).open(encoding='utf8',newline='') as f:return list(csv.DictReader(f))
    spectra=table('spectra.csv');contrasts=table('geometry_contrasts.csv');predictions=table('predictions.csv')
    figdir=folder/'figures';figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.grid':True,'grid.alpha':.25,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':320})
    titles={'A1':'A1: κⱼ=0,1; χₛ=0,5','A2':'A2: κⱼ=0,1; χₛ=10⁻⁴','A3':'A3: κⱼ=1; χₛ=0,01','A4':'A4: κⱼ=10; точная связь'}
    def save(fig,name):
        fig.tight_layout();fig.savefig(figdir/(name+'.pdf'),bbox_inches='tight')
        fig.savefig(figdir/(name+'.png'),dpi=320,bbox_inches='tight');plt.close(fig)
    def curve(ax,identifier,value,label,**kwargs):
        rr=sorted((r for r in spectra if r['id']==identifier),key=lambda r:float(r['beta_deg']))
        ax.plot([float(r['beta_deg']) for r in rr],[float(r[value]) for r in rr],label=label,**kwargs)
    def axes4(ylabel):
        fig,axes=plt.subplots(2,2,figsize=(11,7.2),sharex=True)
        for ax,regime in zip(axes.flat,titles):
            ax.set(title=titles[regime],xlim=(0,90),xlabel='Угол сопряжения β, град',ylabel=ylabel)
        return fig,list(axes.flat)
    for value,name,ylabel in [('omega1','01_thickness_frequency','Первая круговая частота ω₁'),
                              ('Lambda1','02_thickness_normalized','Частотный параметр Λ₁')]:
        fig,axes=axes4(ylabel)
        for ax,regime in zip(axes,titles):
            for epsilon,color in [(.0025,'#0072B2'),(.005,'#009E73'),(.006,'#D55E00')]:
                identifier='B0_'+regime if epsilon==.005 else f'B1_{regime}_e{epsilon:g}'
                curve(ax,identifier,value,f'ε={epsilon:g}',color=color,lw=1.6)
            ax.legend(fontsize=9)
        save(fig,name)
    fig,axes=axes4('Первая круговая частота ω₁')
    colors=plt.get_cmap('viridis')(np.linspace(.08,.90,5))
    for ax,regime in zip(axes,titles):
        for mu,color in zip([-.1,-.01,0,.01,.1],colors):
            identifier='B0_'+regime if mu==.01 else f'B2_{regime}_m{mu:g}'
            curve(ax,identifier,'omega1',f'μ={mu:g}',color=color,lw=1.5)
        ax.legend(fontsize=8,ncol=2)
    save(fig,'03_length_distribution')
    fig,axes=axes4('Изменение ω₁ при aspect 4 → 2, ×10⁻¹⁰')
    for ax,regime in zip(axes,titles):
        rr=[r for r in contrasts if r['regime']==regime and r['group']=='B3']
        ax.plot([float(r['beta_deg']) for r in rr],[float(r['relative_omega1'])*1e10 for r in rr],'.-',lw=.8,ms=3)
        ax.axhline(0,color='black',lw=.7)
    fig.suptitle('Общие расчётные узлы; очень малые разности требуют учёта численной погрешности',fontsize=11)
    save(fig,'04_aspect_difference')
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),sharex=True)
    for ax,regime in zip(axes,['A2','A4']):
        for name,label,style in [('B0_'+regime,'Исходная геометрия','-'),('B1_'+regime+'_e0.006','Только ε=0,006','--'),
                                  ('B2_'+regime+'_m-0.1','Только μ=−0,1',':'),('B4_'+regime,'ε=0,006; μ=−0,1; aspect=2','-.')]:
            curve(ax,name,'omega1',label,lw=1.7,ls=style)
        ax.set(title=titles[regime],xlim=(0,90),xlabel='Угол β, град',ylabel='Первая круговая частота ω₁');ax.legend(fontsize=8)
    save(fig,'05_combined_geometry')
    fig,axes=plt.subplots(2,2,figsize=(11,7.2))
    for ax,regime in zip(axes.flat,titles):
        rr=[r for r in predictions if r['regime']==regime and r['normalization']=='omega90']
        ids=list(dict.fromkeys(r['id'] for r in rr));x=np.arange(len(ids))
        for q,color,marker in [(.9,'#0072B2','o'),(.95,'#D55E00','s')]:
            for i,identifier in enumerate(ids):
                r=next(r for r in rr if r['id']==identifier and float(r['q'])==q)
                if r['combined_error_percent']:
                    ax.plot(i,float(r['combined_error_percent']),marker,color=color,ms=5)
                else:ax.annotate('нет',xy=(i,0),xytext=(0,4),textcoords='offset points',ha='center',color=color,fontsize=8)
            ax.plot([],[],marker,color=color,label=f'{int(q*100)}%')
        ax.axhline(0,color='black',lw=.6)
        labels=[name.replace('_'+regime,'').replace('B1_e','ε=').replace('B2_m','μ=').replace('B3_a','a=').replace('B0','База').replace('B4','Совм.') for name in ids]
        ax.set_xticks(x,labels,rotation=45,ha='right');ax.set(title=titles[regime],ylabel='Ошибка комбинированного угла, %');ax.legend(fontsize=9,loc='upper left')
    save(fig,'06_prediction_errors')
    sources=['spectra.csv','summary.csv','thresholds.csv','geometry.csv','geometry_contrasts.csv','predictions.csv']
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    for name in sources:(figdir/name).write_bytes((folder/name).read_bytes())
    manifest=dict(calculation_mode='plot_only',new_spectral_solves=0,dpi=320,
        sources={name:digest(folder/name) for name in sources},artifacts={p.name:digest(p) for p in figdir.iterdir() if p.is_file()})
    (folder/'figures_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf8')
    print('plot_only',figdir,'zero spectral solves')
