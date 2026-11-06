"""HS-A pure plot_only rendering from accepted CSV; no physics imports."""
from pathlib import Path
import csv
import hashlib
import json
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def render(folder):
    folder=Path(folder); start=time.perf_counter()
    with (folder/'spectra.csv').open(encoding='utf8',newline='') as f:
        spectra=list(csv.DictReader(f))
    with (folder/'summary.csv').open(encoding='utf8',newline='') as f:
        summaries=list(csv.DictReader(f))
    figdir=folder/'figures';figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.grid':True,
                         'grid.alpha':.25,'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'savefig.dpi':320})
    chis=[.0001,.001,.01,.1,.5,'RIGID']
    labels=[r'$\chi_s=10^{-4}$',r'$10^{-3}$',r'$10^{-2}$',r'$0.1$',r'$0.5$','Точная связь']
    colors=plt.get_cmap('tab10').colors[:6]
    def save(fig,name):
        fig.savefig(figdir/(name+'.pdf'),bbox_inches='tight')
        fig.savefig(figdir/(name+'.png'),dpi=320,bbox_inches='tight')
        plt.close(fig)
    for zoom,name in [(False,'01_full_frequency'),(True,'02_small_angles')]:
        fig,axes=plt.subplots(1,3,figsize=(13,4.1),sharey=True)
        for ax,j in zip(axes,[.1,1,10]):
            for chi,label,color in zip(chis,labels,colors):
                rr=[r for r in spectra if float(r['kappa_j'])==j and (r['chi_s']=='RIGID' if chi=='RIGID' else r['chi_s']!='RIGID' and float(r['chi_s'])==chi)]
                rr.sort(key=lambda r:float(r['beta_deg']))
                ax.plot([float(r['beta_deg']) for r in rr],[float(r['omega1']) for r in rr],color=color,label=label,lw=1.6,
                        marker='.',markersize=2.4 if zoom else 1.8)
            ax.set(xlim=(0,6 if zoom else 90),xlabel='Угол сопряжения β, град',title=rf'$\kappa_j={j:g}$')
        axes[0].set_ylabel('Первая круговая частота ω₁')
        axes[0].set_ylim(bottom=0)
        handles,lab=axes[0].get_legend_handles_labels()
        fig.legend(handles,lab,loc='lower center',ncol=6,bbox_to_anchor=(.5,-.035))
        fig.subplots_adjust(bottom=.20,wspace=.12)
        save(fig,name)
    no_drop=all(float(r['drop_to_90_percent'])==0 for r in summaries)
    fig,axes=plt.subplots(1 if no_drop else 2,3,figsize=(13,4.3 if no_drop else 7),squeeze=False,sharex=True)
    for col,j in enumerate([.1,1,10]):
        rr=[r for r in summaries if float(r['kappa_j'])==j]
        rr.sort(key=lambda r:float('inf') if r['chi_s']=='RIGID' else float(r['chi_s']))
        x=np.arange(6)
        axes[0,col].plot(x,[float(r['omega_grid_max']) for r in rr],'o',label='Сеточный максимум',ms=7)
        axes[0,col].plot(x,[float(r['omega90']) for r in rr],'x',label='Частота при 90°',ms=7)
        if not no_drop:
            axes[1,col].bar(x,[float(r['drop_to_90_percent']) for r in rr],color=colors)
            for k,r in enumerate(rr):
                axes[1,col].annotate(f"{float(r['beta_max_grid']):g}°",(k,float(r['drop_to_90_percent'])),
                                     xytext=(0,4),textcoords='offset points',ha='center',fontsize=8)
            axes[1,col].set_ylim(0,max(.012,1.25*max(float(r['drop_to_90_percent']) for r in rr)))
        axes[0,col].set_title(rf'$\kappa_j={j:g}$')
        axes[-1,col].set_xticks(x,['10⁻⁴','10⁻³','10⁻²','0.1','0.5','Точная'])
        axes[-1,col].set_xlabel('Жёсткость χₛ; точная горизонтальная связь')
    axes[0,0].set_ylabel('Круговая частота ω₁')
    if not no_drop:axes[1,0].set_ylabel('Снижение от максимума к 90°, %')
    axes[0,0].legend(fontsize=9)
    fig.suptitle('Для всех 18 конструкций: сеточный максимум при 90°, снижение к 90° равно нулю' if no_drop else
                 'Частотный уровень и снижение к 90°; над столбцами — угол сеточного максимума',fontsize=12)
    fig.tight_layout()
    save(fig,'03_grid_maximum_and_90')
    # Numerical plotting rows retained verbatim, with hashes of source tables.
    for name in ['spectra.csv','summary.csv','thresholds.csv']:
        (figdir/name).write_bytes((folder/name).read_bytes())
    manifest={'calculation_mode':'plot_only','new_eigenvalue_solves':0,'dpi':320,
              'seconds':time.perf_counter()-start,
              'sources':{n:hashlib.sha256((folder/n).read_bytes()).hexdigest() for n in ['spectra.csv','summary.csv','thresholds.csv']},
              'artifacts':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in figdir.iterdir() if p.is_file()}}
    (folder/'figures_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
    print('plot_only',folder,'zero spectral solves')
