#!/usr/bin/env python3
"""Angkor rice aDNA interactive PCA dashboard (plotly, single-file offline HTML).

Reads a smartpca .evec/.eval pair plus an enriched per-sample metadata TSV and
renders one interactive HTML with multiple dropdown views plus static PNG/PDF.

Visual encoding
---------------
- Modern reference: colored BY POPULATION (label from .evec), low opacity,
  small markers, acting as a background coordinate reference.
- Ancient samples:
    * SHAPE = core (CAM2509=circle, CAM23-13=square, CAM23-11=diamond,
      CAM22-08=triangle-up, CAM2201=x).
    * COLOR = continuous Age CE (Viridis) for dated cores; undated cores
      (CAM22-08 / CAM2201) get their own solid colour.
    * Hover: sample_id / core / site / depth / age / libtype / prep / besthit.

Views (top horizontal dropdowns)
--------------------------------
- Panoramic (default): modern + all ancient age/undated points; trajectory and
  QC traces hidden (legendonly).
- Per-core isolation: pick one core.
- QC benchmark: only the multi-library test samples prominent, others faded.
- Trajectory: per-core 100-yr centroid polylines (toggle on).

Static exports (matplotlib)
---------------------------
- {prefix}_panoramic.png, {prefix}_percore.png, {prefix}_pca1_age.png
"""
import argparse
import collections
import sys


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--evec", required=True)
    ap.add_argument("--eval", required=True)
    ap.add_argument("--meta", required=True,
                    help="enriched meta TSV: sample_id base_robot core site "
                         "depth_cm age_CE libtype besthit prep")
    ap.add_argument("--nmarkers", type=int, required=True)
    ap.add_argument("--title", default="angkor full6M (q25)")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--age-min", type=float, default=970.0)
    ap.add_argument("--age-max", type=float, default=2021.0)
    ap.add_argument("--no-png", action="store_true", help="skip matplotlib PNG exports")
    return ap.parse_args()


def load_eval(path):
    vals = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                vals.append(float(line))
    return vals


def load_evec(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.lstrip().startswith("#"):
                continue
            f = line.split()
            if len(f) < 3:
                continue
            rows.append((f[0], [float(x) for x in f[1:-1]], f[-1]))
    return rows


def load_meta(path):
    meta = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            meta[f[0]] = dict(zip(header, f))
    return meta


CORE_SYMBOL = {
    "CAM2509": "circle", "CAM23-13": "square", "CAM23-11": "diamond",
    "CAM22-08": "triangle-up", "CAM2201": "x",
}
CORE_ORDER = ["CAM2509", "CAM23-13", "CAM23-11", "CAM22-08", "CAM2201"]
# solid colours for undated cores (avoid clashing with modern bg)
UNDATED_COLOR = {"CAM22-08": "#e07b39", "CAM2201": "#c2185b"}


def collect(args):
    evals = load_eval(args.eval)
    total = sum(evals) or 1.0
    meta = load_meta(args.meta)
    rows = load_evec(args.evec)
    modern = [(r[0], r[1], r[2]) for r in rows if r[2] != "Ancient"]
    anc_rows = [r for r in rows if r[2] == "Ancient"]

    ancient = collections.defaultdict(list)
    for iid, vals, _ in anc_rows:
        m = meta.get(iid, {})
        core = m.get("core", "?")
        age = None
        if m.get("age_CE"):
            try:
                age = float(m["age_CE"])
            except ValueError:
                age = None
        hover = (f"ID: {iid}<br>core: {core}<br>site: {m.get('site','')}"
                 f"<br>depth: {m.get('depth_cm','')}<br>age: {m.get('age_CE','')}"
                 f"<br>libtype: {m.get('libtype','')}<br>prep: {m.get('prep','')}"
                 f"<br>besthit: {m.get('besthit','')}")
        rec = {"id": iid, "x": vals[0], "y": vals[1], "pc": vals,
               "age": age, "hover": hover, "base": m.get("base_robot", iid),
               "core": core, "libtype": m.get("libtype", "")}
        ancient[core].append(rec)
    cores = [c for c in CORE_ORDER if c in ancient]
    return evals, total, modern, ancient, cores


def robust_range(vals):
    import numpy as np
    a = np.asarray(vals, dtype=float)
    lo, hi = np.nanpercentile(a, 2.0), np.nanpercentile(a, 98.0)
    pad = (hi - lo) * 0.08
    if not np.isfinite(pad) or pad <= 0:
        pad = 0.05
    return float(lo - pad), float(hi + pad)


def trajectory(ancient, core, bin=100):
    pts = [p for p in ancient[core] if p["age"] is not None]
    bins = collections.defaultdict(list)
    for p in pts:
        bins[int(p["age"] // bin) * bin].append(p)
    out = []
    for b in sorted(bins):
        xs = [p["x"] for p in bins[b]]
        ys = [p["y"] for p in bins[b]]
        out.append((f"{b}s", sum(xs) / len(xs), sum(ys) / len(ys)))
    return out


def build_html(args, evals, total, modern, ancient, cores):
    import plotly.graph_objects as go

    # ---- modern traces: one per population (coloured), low opacity ----
    pop_traces = []
    pop_list = sorted({lab for _, _, lab in modern})
    palette = ["#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
               "#ffd8a2", "#98d8c8", "#f7b6d2", "#c7c7c7", "#9edae5",
               "#dbdb8d", "#ff9f40", "#bcbd22", "#17becf", "#8c564b",
               "#e377c2", "#7f7f7f", "#bcbd22", "#d62728", "#2ca02c"]
    for i, lab in enumerate(pop_list):
        xs = [r[1][0] for r in modern if r[2] == lab]
        ys = [r[1][1] for r in modern if r[2] == lab]
        pop_traces.append(go.Scatter(
            x=xs, y=ys, mode="markers", name=lab,
            marker=dict(color=palette[i % len(palette)], size=5, opacity=0.30,
                        line=dict(width=0)), hoverinfo="skip", showlegend=True))
    n_pop = len(pop_traces)

    # ---- ancient traces: per core (age-coloured) + undated solid ----
    ancient_traces = []
    core_trace_idx = {}   # core -> list of trace indices for that core's age pts
    undated_idx = {}      # core -> trace index (solid colour)
    for core in cores:
        pts = ancient[core]
        aged = [p for p in pts if p["age"] is not None]
        noage = [p for p in pts if p["age"] is None]
        sym = CORE_SYMBOL.get(core, "circle")
        if aged:
            core_trace_idx[core] = len(ancient_traces)
            ancient_traces.append(go.Scatter(
                x=[p["x"] for p in aged], y=[p["y"] for p in aged],
                mode="markers", name=core,
                marker=dict(size=9, symbol=sym,
                            color=[p["age"] for p in aged],
                            colorscale="Viridis", cmin=args.age_min, cmax=args.age_max,
                            colorbar=dict(title="Age CE", x=1.02, xanchor="left"),
                            line=dict(width=0.6, color="white")),
                customdata=[p["hover"] for p in aged],
                hovertemplate="%{customdata}<extra></extra>", showlegend=True))
        else:
            core_trace_idx[core] = None
        if noage:
            col = UNDATED_COLOR.get(core, "#888888")
            undated_idx[core] = len(ancient_traces)
            ancient_traces.append(go.Scatter(
                x=[p["x"] for p in noage], y=[p["y"] for p in noage],
                mode="markers", name=f"{core} (undated)",
                marker=dict(size=9, symbol=sym, color=col,
                            line=dict(width=0.6, color="white")),
                customdata=[p["hover"] for p in noage],
                hovertemplate="%{customdata}<extra></extra>", showlegend=True))
    n_anc = len(ancient_traces)

    # ---- trajectory traces (legendonly by default) ----
    traj_traces = []
    for core in cores:
        tr = trajectory(ancient, core)
        if len(tr) < 2:
            continue
        traj_traces.append(go.Scatter(
            x=[p[1] for p in tr], y=[p[2] for p in tr],
            mode="lines+markers", name=f"{core} traj",
            line=dict(color="black", width=2), marker=dict(size=8, color="red"),
            hoverinfo="skip", visible="legendonly", showlegend=True))
    n_traj = len(traj_traces)

    # ---- QC markers + lines (legendonly by default) ----
    libgroup = collections.defaultdict(list)
    for core in cores:
        for p in ancient[core]:
            if p["libtype"] in ("SG", "C1", "C2", "pooled_nobesthit", "pooled_besthit"):
                libgroup[p["base"]].append(p)
    order = {"SG": 0, "C1": 1, "C2": 2, "pooled_nobesthit": 3, "pooled_besthit": 4}
    qc_marker_traces = []
    qc_line_traces = []
    qc_bases = []
    for base, pts in libgroup.items():
        if len(pts) < 2:
            continue
        pts = sorted(pts, key=lambda p: order.get(p["libtype"], 9))
        qc_bases.append(base)
        qc_marker_traces.append(go.Scatter(
            x=[p["x"] for p in pts], y=[p["y"] for p in pts],
            mode="markers", name=base, visible="legendonly",
            marker=dict(size=11, color="black", symbol="circle",
                        line=dict(color="white", width=1)),
            customdata=[p["hover"] for p in pts],
            hovertemplate="%{customdata}<extra></extra>", showlegend=False))
        qc_line_traces.append(go.Scatter(
            x=[p["x"] for p in pts] + [pts[0]["x"]],
            y=[p["y"] for p in pts] + [pts[0]["y"]],
            mode="lines", name=f"{base} QC line", visible="legendonly",
            line=dict(color="rgba(0,0,0,0.6)", width=1.5, dash="dot"),
            hoverinfo="skip", showlegend=False))
    n_qc_m = len(qc_marker_traces)
    n_qc_l = len(qc_line_traces)

    # combined trace layout:
    # [pop_traces (n_pop)] [ancient_traces (n_anc)] [traj (n_traj)] [qc_m (n_qc_m)] [qc_l (n_qc_l)]
    pop_vis = [True] * n_pop
    anc_all = [1] * n_anc

    def vis(pop_on, anc, traj_on, qc_m_on, qc_l_on):
        v = ([True] * n_pop if pop_on else [False] * n_pop)
        v += ([1] * n_anc if anc else [0] * n_anc)
        v += ([1] * n_traj if traj_on else [0] * n_traj)
        v += ([1] * n_qc_m if qc_m_on else [0] * n_qc_m)
        v += ([1] * n_qc_l if qc_l_on else [0] * n_qc_l)
        return v

    # ---- view buttons ----
    # panoramic
    b_pano = dict(label="Panoramic", method="update",
                  args=[{"visible": vis(True, True, False, False, False)}])
    # per-core isolation
    core_buttons = [dict(label="All cores", method="update",
                         args=[{"visible": vis(True, True, False, False, False)}])]
    for core in cores:
        mask = [0] * n_anc
        ci = core_trace_idx.get(core)
        ui = undated_idx.get(core)
        if ci is not None:
            mask[ci] = 1
        if ui is not None:
            mask[ui] = 1
        core_buttons.append(dict(label=f"Core: {core}", method="update",
                                 args=[{"visible": vis(True, mask, False, False, False)}]))
    # QC benchmark: modern + qc markers + qc lines, ancient dimmed
    b_qc = dict(label="QC benchmark", method="update",
                args=[{"visible": vis(True, False, False, True, True)}])
    # trajectory
    b_traj = dict(label="Trajectory", method="update",
                  args=[{"visible": vis(True, True, True, False, False)}])

    # ---- robust axis range ----
    all_x = [p[1][0] for p in modern] + [p["x"] for c in cores for p in ancient[c]]
    all_y = [p[1][1] for p in modern] + [p["y"] for c in cores for p in ancient[c]]
    xr = robust_range(all_x)
    yr = robust_range(all_y)

    all_traces = pop_traces + ancient_traces + traj_traces + qc_marker_traces + qc_line_traces
    fig = go.Figure(data=all_traces)
    fig.update_layout(
        title=dict(text=f"{args.title} | markers={args.nmarkers}", font=dict(size=15)),
        hovermode="closest",
        xaxis=dict(title=f"PC1 ({evals[0]/total*100:.2f}%)", range=xr),
        yaxis=dict(title=f"PC2 ({evals[1]/total*100:.2f}%)", range=yr),
        height=720, width=1100,
        margin=dict(l=60, r=70, t=90, b=60),
        legend=dict(orientation="h", y=1.02, x=0, font=dict(size=10)),
        updatemenus=[
            dict(buttons=[b_pano, b_qc, b_traj] + core_buttons,
                 direction="down", showactive=True,
                 x=0.0, y=1.10, xanchor="left", yanchor="top",
                 bgcolor="rgba(240,240,240,0.9)"),
        ],
    )

    html_out = f"{args.out_prefix}.dashboard.html"
    fig.write_html(html_out, include_plotlyjs="cdn", full_html=True)
    print(f"wrote {html_out}")


def build_pngs(args, evals, total, modern, ancient, cores):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pops = sorted({r[2] for r in modern})
    pop_color = {p: plt.get_cmap("tab20")(i % 20) for i, p in enumerate(pops)}

    def draw(ax, data):
        for p in pops:
            xs = [r[1][0] for r in modern if r[2] == p]
            ys = [r[1][1] for r in modern if r[2] == p]
            ax.scatter(xs, ys, s=8, color=pop_color[p], alpha=0.30, linewidths=0)
        for c, pts in data.items():
            aged = [p for p in pts if p["age"] is not None]
            noage = [p for p in pts if p["age"] is None]
            if aged:
                ax.scatter([p["x"] for p in aged], [p["y"] for p in aged],
                           c=[p["age"] for p in aged], cmap="viridis", s=40,
                           vmin=args.age_min, vmax=args.age_max,
                           edgecolor="white", linewidth=0.5)
            if noage:
                col = UNDATED_COLOR.get(c, "#888888")
                ax.scatter([p["x"] for p in noage], [p["y"] for p in noage],
                           c=col, s=40, edgecolor="white", linewidth=0.5)
        ax.set_xlabel(f"PC1 ({evals[0]/total*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({evals[1]/total*100:.1f}%)")

    fig, ax = plt.subplots(figsize=(9, 7))
    draw(ax, dict(ancient))
    sm = plt.cm.ScalarMappable(cmap="viridis",
                               norm=plt.Normalize(args.age_min, args.age_max))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Age CE")
    ax.set_title(f"{args.title} | markers={args.nmarkers}")
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_panoramic.png", dpi=150)
    print(f"wrote {args.out_prefix}_panoramic.png")

    fig, axes = plt.subplots(1, len(cores), figsize=(6 * len(cores), 5.5))
    if len(cores) == 1:
        axes = [axes]
    for ax, core in zip(axes, cores):
        draw(ax, {core: ancient.get(core, [])})
        ax.set_title(core)
    fig.suptitle(f"{args.title} | markers={args.nmarkers}")
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_percore.png", dpi=150)
    print(f"wrote {args.out_prefix}_percore.png")

    fig, ax = plt.subplots(figsize=(9, 6))
    for core in cores:
        pts = [(p["age"], p["pc"][0]) for p in ancient.get(core, [])
               if p["age"] is not None]
        if not pts:
            continue
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts], marker="o",
                label=core, linewidth=1.5)
    ax.set_xlabel("Age CE")
    ax.set_ylabel(f"PC1 ({evals[0]/total*100:.1f}%)")
    ax.set_title(f"{args.title} | PC1 vs Age")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_pca1_age.png", dpi=150)
    print(f"wrote {args.out_prefix}_pca1_age.png")


def main():
    args = parse_args()
    evals, total, modern, ancient, cores = collect(args)
    build_html(args, evals, total, modern, ancient, cores)
    if not args.no_png:
        try:
            build_pngs(args, evals, total, modern, ancient, cores)
        except ImportError as e:
            print(f"WARNING: matplotlib not available, skipping PNG ({e})",
                  file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
