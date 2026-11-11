"""Publication figures from immutable RL/RR BETA tables; no spectral evaluations.

Separate presentation I/O preserves the source-hashed calculation entry point
and its diagnostic defaults. Reuses map_rl_full_rr.saved_table for verified CSV
reading. Windows come from the configuration, never from interpolation.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.analysis import map_rl_full_rr as archive
from scripts.lib import frequency_geometry_style as line_geometry

COLORS = ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
          "#9467bd", "#8c564b", "#b85d93", "#666666")


def read_maps(folder, config):
    manifest, rows = archive.saved_table(folder)
    if manifest["status"] != "PASS":
        raise ValueError("Unresolved bundle: publication output stopped")
    specs = manifest["identity"]["maps"]
    k_plot, k_guard = archive.prefix_ranks(config)
    if archive.prefix_ranks(manifest["identity"]) != (k_plot,k_guard):
        raise ValueError("Configuration prefix does not describe the saved bundle")
    if specs != config["maps"] or manifest["identity"]["material"] != config["material"]:
        raise ValueError("Configuration does not describe the saved bundle")
    if any(s["mode"] != "BETA" for s in specs):
        raise ValueError("This presentation entry point accepts BETA maps only")
    if any(p["spectrum_semantics"] != "sorted_positions"
           for p in manifest["identity"]["local_policies"].values()):
        raise ValueError("Sorted-position tables required")
    indexed = {}
    for row in rows:
        key = (row["map"], float(row["x"]), int(row["rank"]))
        if key in indexed:
            raise ValueError(f"Duplicate CSV entry: {key}")
        if row["status"] != "PASS" or not math.isfinite(float(row["Lambda"])):
            raise ValueError(f"Unresolved CSV entry: {key}")
        indexed[key] = float(row["Lambda"])
    result = []
    for spec in specs:
        x = spec["grid"]
        window = config["grid_provenance"]["local_windows"][spec["name"]]
        if not x[0] <= window[0] < window[1] <= x[-1]:
            raise ValueError("Invalid local display window")
        y = [[indexed[(spec["name"], b, rank)] for rank in range(1, k_guard+1)] for b in x]
        if any(any(a > b for a, b in zip(v, v[1:])) for v in y):
            raise ValueError("CSV frequencies are not sorted")
        result.append(dict(spec=spec, x=x, y=y, window=window, material=config["material"], K_plot=k_plot))
    return manifest, result


def sampled_gaps(maps, config):
    from decimal import Decimal
    step = Decimal(str(config["grid_provenance"]["analysis_step"]))
    local_step = Decimal(str(config["grid_provenance"]["local_step"]))
    rows = []
    for data in maps:
        i = min(range(len(data["x"])), key=lambda i: data["y"][i][1]-data["y"][i][0])
        beta, y = data["x"][i], data["y"][i]
        low, high = data["window"]
        b = Decimal(str(beta))
        main = (b-Decimal(str(data["x"][0]))) % step == 0
        local = low <= beta <= high and (b-Decimal(str(low))) % local_step == 0
        rows.append(dict(map=data["spec"]["name"], mu=data["spec"]["fixed"]["mu"],
            beta_deg=beta, Lambda_1=y[0], Lambda_2=y[1], Delta_Lambda_12=y[1]-y[0],
            grid="main+local" if main and local else "main" if main else "local",
            interior_to_local_window=low < beta < high))
    return rows


def draw(ax, data, ranks, *, zoom=False, ymax=None):
    selected = [i for i, x in enumerate(data["x"])
                if not zoom or data["window"][0] <= x <= data["window"][1]]
    x = [data["x"][i] for i in selected]
    geometry = line_geometry.FrequencyGeometry.from_spec(data["spec"], data["material"], x)
    for rank in ranks:
        line_geometry.plot_frequency_curve(ax, geometry, [data["y"][i][rank-1] for i in selected],
                color=COLORS[rank-1], linewidth=1.35, marker=None)
    ax.set_xlabel(r"$\beta$, градусы", labelpad=3)
    ax.set_ylabel(r"$\Lambda$", rotation=0, labelpad=13)
    ax.set_xlim(data["window"] if zoom else (data["x"][0], data["x"][-1]))
    if ymax is not None:
        ax.set_ylim(0, ymax)
    else:
        ax.margins(y=.10)
    ax.grid(False, which="both")
    ax.tick_params(direction="in", top=True, right=True, length=3, width=.7)
    for spine in ax.spines.values():
        spine.set_linewidth(.7)
    # Original observations stay unchanged; any extra vertex is display-only.
    assert not any(line.get_visible() for line in ax.get_xgridlines()+ax.get_ygridlines())
    return dict(map=data["spec"]["name"], ranks=list(ranks), nodes=x,
                no_markers=True, no_grid=True, straight_segments=True, geometry=geometry.diagnostic())


def render(folder, config_path, *, dpi=260, font_size=12):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    config = archive.read(config_path)
    manifest, maps = read_maps(folder, config)
    k_plot, _ = archive.prefix_ranks(config)
    count_word = "шесть" if k_plot==6 else "восемь"
    font = font_manager.findfont("Times New Roman", fallback_to_default=False)
    identity = dict(source_bundle=folder.resolve().relative_to(ROOT).as_posix(),
        source_manifest_sha256=archive.sha(folder/"manifest.json"),
        source_csv_sha256=archive.sha(folder/"spectrum.csv"),
        config=config_path.resolve().relative_to(ROOT).as_posix(),
        config_sha256=archive.sha(config_path),
        renderer_sha256=archive.sha(__file__), reader_sha256=archive.sha(archive.__file__),
        geometry_sources=line_geometry.provenance(),
        matplotlib=matplotlib.__version__, font=font, font_sha256=archive.sha(font),
        dpi=dpi, font_size=font_size, colors=COLORS, spectrum_semantics="sorted_positions",
        calculation_mode="plot_only", root_calculations=0)
    output = ROOT/"results/rl_full_rr_maps/figures"/(folder.name+"_publication_"+archive.digest(identity)[:12])
    if (output/"manifest.json").exists():
        archive.verify_manifest(output)
        return output
    if output.exists():
        raise FileExistsError("Incomplete presentation folder; preserve it and use changed style settings")
    output.mkdir(parents=True)
    geometry_diagnostics = line_geometry.write_geometry_diagnostics(output, {
        d["spec"]["name"]: line_geometry.FrequencyGeometry.from_spec(d["spec"], d["material"])
        for d in maps})
    gaps = sampled_gaps(maps, config)
    with (output/"sampled_gaps.csv").open("x", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(gaps[0])); writer.writeheader(); writer.writerows(gaps)
    ymax = math.ceil(max(max(y[:k_plot]) for d in maps for y in d["y"]))
    plots = []
    epsilon = config["material"]["epsilon"]
    def params(d, letter=None):
        prefix = f"({letter})  " if letter else ""
        return prefix+rf"$\mu={d['spec']['fixed']['mu']:g},\quad\tau={d['spec']['fixed']['tau']:g},\quad\varepsilon={epsilon:g}$"
    def save(fig, name, panels, kind):
        caption = line_geometry.add_geometry_caption(fig)
        for ext in ("png", "pdf"):
            fig.savefig(output/(name+"."+ext), dpi=dpi, facecolor="white", bbox_inches="tight")
        plots.append(dict(name=name, kind=kind, panels=panels, geometry_caption=caption))
        plt.close(fig)
    with plt.rc_context({"font.family":"serif", "font.serif":["Times New Roman"],
            "font.size":font_size, "mathtext.fontset":"stix", "pdf.fonttype":42,
            "axes.grid":False, "figure.facecolor":"white", "axes.facecolor":"white",
            "path.simplify":False}):
        for d in maps:
            for zoom in (False, True):
                fig, ax = plt.subplots(figsize=(6.3, 4.9))
                fig.subplots_adjust(left=.12, right=.97, top=.96, bottom=.27)
                panel = draw(ax, d, (1,2) if zoom else range(1,k_plot+1), zoom=zoom, ymax=None if zoom else ymax)
                fig.text(.55, .145, params(d), ha="center")
                caption = ("Первые две упорядоченные безразмерные собственные частоты.\n"
                           "Увеличенный фрагмент угловой зависимости." if zoom else
                           f"Первые {count_word} упорядоченных безразмерных собственных частот\n"
                           "в зависимости от угла сопряжения.")
                fig.text(.5, .038, caption, ha="center", va="bottom", fontsize=font_size-1)
                save(fig, d["spec"]["name"]+("_pair" if zoom else "_all"), [panel], "pair" if zoom else "all")
        columns = min(3, len(maps)); rows = math.ceil(len(maps)/columns)
        fig, axes = plt.subplots(rows, columns, figsize=(columns*4.2, rows*3.25+1), squeeze=False)
        fig.subplots_adjust(left=.055, right=.98, top=.97, bottom=.14, hspace=.62, wspace=.27)
        panels=[]
        for i, ax in enumerate(axes.flat):
            if i >= len(maps):
                ax.set_visible(False); continue
            d=maps[i]; panels.append(draw(ax,d,range(1,k_plot+1),ymax=ymax))
            ax.text(.5,-.36,params(d,chr(ord('a')+i)),transform=ax.transAxes,ha="center",fontsize=font_size-1)
        fig.text(.5,.04,f"Первые {count_word} упорядоченных безразмерных собственных частот в зависимости от угла сопряжения.\n"
                 "Обзор случаев с различными значениями параметра неравенства длин.",ha="center",fontsize=font_size)
        save(fig,"overview",panels,"overview")
    gallery = ["# Угловые зависимости: публикационное оформление", "",
        "Линии соединяют рассчитанные узлы без сглаживания. Цвета обозначают",
        f"упорядоченные позиции 1–{k_plot} снизу вверх; продолжение форм не выполнялось.",
        "Обзорный лист служит для сопоставления случаев, а не окончательным рисунком главы.", "",
        f"Исходные [CSV](../../{folder.name}/spectrum.csv) и [manifest](../../{folder.name}/manifest.json).",
        "[Минимальные зазоры на рассчитанных узлах](sampled_gaps.csv).", ""]
    for p in plots:
        gallery.extend([f"## {p['name']}", "", f"[PNG]({p['name']}.png) · [PDF]({p['name']}.pdf)", "",
                        f"![{p['name']}]({p['name']}.png)", ""])
    (output/"README.md").write_text("\n".join(gallery),encoding="utf-8")
    artifacts={p.name:archive.sha(p) for p in output.iterdir() if p.is_file()}
    archive.publish(output/"manifest.json",dict(status="PASS",identity=identity,plots=plots,
        geometry_diagnostics=geometry_diagnostics,
        artifacts=artifacts, shared_full_ylim=[0,ymax], pair_limits="local windows, individually scaled Lambda",
        source_statistics=manifest["statistics"]))
    return output


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle",type=Path,required=True)
    parser.add_argument("--config",type=Path,required=True)
    parser.add_argument("--dpi",type=int,default=260)
    parser.add_argument("--font-size",type=float,default=12)
    args=parser.parse_args()
    print(render(args.bundle,args.config,dpi=args.dpi,font_size=args.font_size))


if __name__ == "__main__":
    main()
