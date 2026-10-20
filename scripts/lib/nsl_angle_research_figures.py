"""Exactly three ANGLE-C figures from prepared tables; no model/solver import.

Fixed research figure helper, not a general plotting application. Figure data
are prepared from accepted saved modes/roots by the diagnostic entry point.
"""
from pathlib import Path
import csv
import hashlib
import json
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


TABLES = ("figure1_maxima.csv", "figure2_rotations.csv", "figure2_pairs.csv",
          "figure3_envelopes.csv", "figure3_modes.csv")
NAMES = ("figure1_angle_maxima", "figure2_joint_rotations", "figure3_modal_envelopes")
COLORS = ("#2166ac", "#d95f02")


def read_rows(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def render(data, destination, dpi=220):
    """Read CSV only; never evaluate or recover an eigenfunction."""
    data, destination = Path(data), Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    with plt.rc_context({"font.family":"DejaVu Sans", "font.size":11,
                         "axes.spines.top":False, "axes.spines.right":False,
                         "pdf.fonttype":42}):
        rows = read_rows(data/TABLES[0])
        fig, ax = plt.subplots(figsize=(8.2, 5.2), layout="constrained")
        for case, color, marker, offset in zip(("MIRROR","SAME_LOCAL"), COLORS, ("o","s"), (-24,24)):
            selected = [r for r in rows if r["case"] == case]
            x = [float(r["beta_deg"]) for r in selected]
            y = [float(r["zeta_pred_max"])*1000 for r in selected]
            ax.scatter(x, y, color=color, marker=marker, s=72, label=case, zorder=3)
            for r, xx, yy in zip(selected, x, y):
                ax.annotate(f'k={r["position"]}', (xx,yy), xytext=(offset,12 if case=="SAME_LOCAL" else -17),
                            textcoords="offset points", ha="center", color=color, fontsize=10)
        ax.set(xticks=[15,30,45,60,75], xlim=(9,81), ylim=(.75,2.16),
               xlabel=r"Угол сопряжения $\beta$, градусы",
               ylabel=r"$10^3\,\zeta_{\max,\mathrm{pred}}$",
               title="Максимум упругого прогноза среди позиций 1–6"
                     "\n"+r"$\kappa_\theta=1,\ d_\theta=0.001$")
        ax.grid(alpha=.2); ax.legend(loc="upper right")
        fig.text(.51,.015,"Отдельные сохранённые углы; k — отсортированная позиция, не номер ветви",
                 ha="center",fontsize=9)
        fig.get_layout_engine().set(h_pad=.16, w_pad=.08, rect=(0,.055,1,.945))
        save(fig, destination, NAMES[0], dpi)

        rotations = read_rows(data/TABLES[1]); pairs = read_rows(data/TABLES[2])
        fig, axes = plt.subplots(2,2,figsize=(11.2,7.3), sharex=True,
                                 sharey="row", layout="constrained")
        for i, position in enumerate((4,6)):
            group = [r for r in rotations if int(r["position"]) == position]
            peak = max(abs(float(r["psi"])) for r in group)
            for j, case in enumerate(("MIRROR","SAME_LOCAL")):
                ax = axes[i,j]
                selected = [r for r in group if r["case"] == case]
                for arm, color, style, marker in zip((1,2),COLORS,("-","--"),("o","x")):
                    values = [r for r in selected if int(r["arm"]) == arm]
                    ax.plot([float(r["xi"]) for r in values],[float(r["psi"]) for r in values],
                            color=color,ls=style,lw=1.8,label=f"Плечо {arm}")
                    end = next(r for r in values if float(r["xi"]) == 1.)
                    ax.plot(1.,float(end["psi"]),marker=marker,color=color,ms=8,
                            mfc="none",mew=1.7,clip_on=False,zorder=5)
                pair = next(r for r in pairs if int(r["position"])==position and r["case"]==case)
                if case=="MIRROR":
                    text = r"$\Delta\psi=0$ точно; $\zeta=0$"
                else:
                    text = (r"$\Delta\psi=$"+f'{float(pair["delta_plot"]):+.5f}'
                            +"\n"+r"$\zeta=$"+f'{float(pair["actual_zeta"]):.6g}')
                ax.text(.03,.97,text,transform=ax.transAxes,va="top",fontsize=10,
                        bbox=dict(facecolor="white",edgecolor="none",alpha=.85))
                ax.set(title=f"{case} · упругая позиция {position}",
                       ylim=(-1.42*peak,1.42*peak),xlim=(0,1.035))
                ax.axhline(0,color=".65",lw=.6);ax.grid(alpha=.16)
                ax.legend(loc="lower left",fontsize=9)
                if j==0:ax.set_ylabel(r"$\psi_i$ при $\mathcal{M}=1$")
                if i==1:ax.set_xlabel(r"$x_i/L_i$  (узел при 1)")
        fig.suptitle(r"Повороты сохранённых упругих форм: $\beta=30^\circ$, $\kappa_\theta=1$"
                     "\nЗнаки согласованы; масштаб одинаков в каждой строке",fontsize=13)
        save(fig,destination,NAMES[1],dpi)

        envelopes=read_rows(data/TABLES[3]); modes=read_rows(data/TABLES[4])
        fig,ax=plt.subplots(figsize=(8.4,5.4),layout="constrained")
        for mode,color in zip(modes,COLORS):
            values=[r for r in envelopes if r["case"]==mode["case"]]
            label=(f'{mode["case"]} 5: '+r"$\zeta=$"+f'{float(mode["zeta"]):.7f}'
                   +", "+r"$N_{1/e}=$"+f'{float(mode["N_1e"]):.2f}')
            ax.plot([float(r["N"]) for r in values],[float(r["amplitude"]) for r in values],
                    color=color,lw=2,label=label)
            ax.plot(float(mode["N_1e"]),np.exp(-1),marker="o",color=color,ms=5)
        ax.axhline(np.exp(-1),ls=":",color=".55",lw=1)
        ax.set(xlim=(0,200),ylim=(0,1.04),xlabel=r"Число собственных периодов $N$",
               ylabel=r"Относительная амплитуда $A/A_0$",
               title=r"Огибающие отдельных собственных мод, $\beta=15^\circ$"
                     "\n"+r"$\kappa_\theta=1,\ d_\theta=0.001$; фактические комплексные корни")
        ax.legend(loc="upper right",fontsize=9);ax.grid(alpha=.18)
        fig.text(.52,.015,"Одинаковое N соответствует разному времени; это не вынужденный отклик",
                 ha="center",fontsize=9)
        fig.get_layout_engine().set(h_pad=.16, w_pad=.08, rect=(0,.055,1,.945))
        save(fig,destination,NAMES[2],dpi)
    manifest=dict(status="PASS",new_eigenvalue_solves=0,figure_count=3,
        matplotlib=matplotlib.__version__,renderer_sha256=digest(__file__),
        table_hashes={n:digest(data/n) for n in TABLES},
        artifact_hashes={n+ext:digest(destination/(n+ext)) for n in NAMES for ext in (".png",".pdf")},
        render_seconds=time.perf_counter()-started)
    (destination/"figure_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    return manifest


def save(fig, destination, name, dpi):
    fig.savefig(destination/(name+".png"),dpi=dpi)
    fig.savefig(destination/(name+".pdf"))
    plt.close(fig)
