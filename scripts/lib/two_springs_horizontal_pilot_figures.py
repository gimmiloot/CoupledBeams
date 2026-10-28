"""Plot-only HS-P1 figures: reads accepted CSV, never imports a solver."""
from pathlib import Path
import csv
import hashlib
import json
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


def read(path):
    with path.open(encoding='utf8',newline='') as f:
        return list(csv.DictReader(f))


def save_csv(path,rows):
    with path.open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)


def render(folder):
    start=time.perf_counter();folder=Path(folder)
    manifest=json.loads((folder/'manifest.json').read_text())
    if manifest['status']!='PASS':raise ValueError('Only accepted complete pilot can be rendered')
    analysis=json.loads((folder/'analysis_manifest.json').read_text())
    for name in ['spectra.csv','energy_distribution.csv','small_angle_comparison.csv','small_angle_coefficients.csv']:
        expected=(manifest if name=='spectra.csv' else analysis)['artifacts'][name]
        if hashlib.sha256((folder/name).read_bytes()).hexdigest()!=expected:
            raise ValueError('Changed plotting input: '+name)
    rows=read(folder/'spectra.csv');energy=read(folder/'energy_distribution.csv')
    asym=read(folder/'small_angle_comparison.csv');coeff=read(folder/'small_angle_coefficients.csv')
    figures=folder/'figures';figures.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':11,
                         'axes.titlesize':11,'legend.fontsize':9,'axes.spines.top':False,
                         'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':320,
                         'axes.grid':True,'grid.alpha':.22})
    colors={'H1':'#0072B2','H2':'#777777','H3':'#009E73','H4':'#E69F00',
            'H5':'#D55E00','H6':'#202020','H7':'#8C56A2'}
    def series(name):
        return sorted([r for r in rows if r['id']==name],key=lambda r:float(r['beta_deg']))
    def line(ax,name,label):
        r=series(name)
        ax.plot([float(x['beta_deg']) for x in r],[float(x['omega1']) for x in r],
                label=label,color=colors[name],ls='--' if name=='H6' else '-',
                marker='o',ms=3,mfc='white' if name=='H6' else colors[name],lw=1.5)
    def axes(ax):
        ax.set_xlabel(r'Угол сопряжения $\beta$, град')
        ax.set_ylabel(r'Первая круговая частота $\omega_1$')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,-3),useMathText=True)
    def save(fig,name,data):
        fig.savefig(figures/(name+'.pdf'),bbox_inches='tight')
        fig.savefig(figures/(name+'.png'),dpi=320,bbox_inches='tight')
        save_csv(figures/(name+'.csv'),data)
        plt.close(fig)

    fig,ax=plt.subplots(figsize=(8.2,4.6),layout='constrained')
    for name,label in [('H2',r'H2: $\kappa_s=0$'),('H3',r'H3: $\kappa_s=1$'),
                       ('H4',r'H4: $\kappa_s=400$'),('H5',r'H5: $\kappa_s=20000$'),
                       ('H6','H6: точная горизонтальная связь')]:line(ax,name,label)
    axes(ax);ax.set_xlim(0,90);ax.set_ylim(0,.0093);ax.set_xticks([0,15,30,45,60,75,90])
    ax.set_title(r'Горизонтальная опора при $\kappa_j=0{,}1$')
    ax.legend(loc='lower right',framealpha=.95)
    save(fig,'01_horizontal_stiffness',[r for r in rows if r['id'] in ['H2','H3','H4','H5','H6']])

    fig,axs=plt.subplots(1,2,figsize=(10,4.1),layout='constrained')
    for ax in axs:
        for name,label in [('H1',r'H1: $\kappa_j=0$'),('H5',r'H5: $\kappa_j=0{,}1$'),('H7',r'H7: $\kappa_j=10$')]:line(ax,name,label)
        axes(ax)
    axs[0].set_xlim(0,5);axs[0].set_ylim(0,.0116);axs[0].set_title('Начальный участок');axs[0].legend(loc='lower right')
    axs[1].set_xlim(5,90);axs[1].set_ylim(.0073,.0116);axs[1].set_title('Участок 5–90°')
    fig.suptitle(r'Внутреннее сопряжение при $\kappa_s=20000$')
    save(fig,'02_joint_stiffness',[r for r in rows if r['id'] in ['H1','H5','H7']])

    fig,axs=plt.subplots(2,3,figsize=(10.4,6.2),layout='constrained')
    plotted=[]
    for ax,name in zip(axs.ravel(),['H1','H3','H4','H5','H6']):
        coef=next(r for r in coeff if r['id']==name)
        c=float(coef['coefficient_beta_rad_squared']);w0=float(coef['omega0'])
        beta=np.linspace(0,.2,101)
        ax.plot(beta,c*np.radians(beta)**2,color='black',ls='--',lw=1.25,label='Малоугловая формула')
        rr=[r for r in asym if r['id']==name and float(r['beta_deg'])<=.2]
        ax.plot([0]+[float(r['beta_deg']) for r in rr],[0]+[float(r['omega2_full'])-w0*w0 for r in rr],
                'o',color=colors[name],ms=4,label='Полная модель')
        ax.set_title(name);ax.set_xlabel(r'$\beta$, град');ax.set_ylabel(r'$\omega_1^2(\beta)-\omega_1^2(0)$')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useMathText=True)
        ax.set_xticks([0,.05,.1,.2]);plotted.extend(rr)
    handles,labels=axs[0,0].get_legend_handles_labels()
    axs[1,2].axis('off');axs[1,2].legend(handles,labels,loc='upper left',frameon=False)
    axs[1,2].text(.03,.55,'Углы в формулах — радианы.\nКоэффициенты не подгонялись.\nУ H1 исходная частота равна нулю.\nПоказан только локальный участок.',
                  transform=axs[1,2].transAxes,ha='left',va='top',fontsize=9,linespacing=1.6)
    fig.suptitle('Локальная проверка малоуглового приближения')
    save(fig,'03_small_angle',plotted)

    fig,axs=plt.subplots(2,1,figsize=(10,7),layout='constrained',sharex=True,
                         gridspec_kw={'height_ratios':[1.5,1]})
    x=np.arange(len(energy));bottom=np.zeros(len(energy))
    components=[('bending_fraction','Изгиб','#0072B2'),('shear_fraction','Сдвиг','#56B4E9'),
                ('axial_fraction','Продольная деформация','#009E73'),
                ('joint_fraction','Внутренняя пружина','#CC79A7'),
                ('horizontal_fraction','Горизонтальная пружина','#E69F00')]
    for field,label,color in components:
        value=np.array([float(r[field]) for r in energy])*100
        axs[0].bar(x,value,bottom=bottom,label=label,color=color,width=.7);bottom+=value
    axs[0].set_ylim(0,100);axs[0].set_ylabel('Доля потенциальной энергии, %')
    axs[0].legend(loc='upper center',bbox_to_anchor=(.5,1.24),ncol=3,frameon=False)
    axs[1].bar(x,[100*float(r['arm1_mass_fraction']) for r in energy],width=.7,color='#555555')
    axs[1].set_ylabel('Вклад первого стержня\nв массовую норму, %');axs[1].set_ylim(0,25)
    axs[1].set_xticks(x,[r['id']+'\n'+f"{float(r['beta_deg']):g}"+'°' for r in energy])
    axs[1].set_xlabel('Конструкция и угол сопряжения; полная массовая норма равна единице')
    save(fig,'04_energy_distribution',energy)

    def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    record={'calculation_mode':'plot_only','new_eigenvalue_solves':0,
            'render_seconds':time.perf_counter()-start,'dpi':320,'labels_language':'Russian',
            'renderer_sha256':digest(__file__ if isinstance(__file__,Path) else Path(__file__)),
            'inputs':{name:digest(folder/name) for name in ['spectra.csv','energy_distribution.csv','small_angle_comparison.csv','small_angle_coefficients.csv']},
            'artifacts':{p.relative_to(folder).as_posix():digest(p) for p in figures.iterdir() if p.is_file()}}
    (folder/'figures_manifest.json').write_text(json.dumps(record,indent=2),encoding='utf8')
    print('plot_only: 4 PDF / 4 PNG / 4 CSV; zero eigenvalue solves',folder)
