#!/usr/bin/env python3
"""Angkor rice aDNA interactive PCA dashboard (plotly, single-file offline HTML).

Reads a smartpca .evec/.eval pair plus an enriched per-sample metadata TSV and
renders one interactive HTML with 4 dropdown views plus static PNG/PDF exports.

Visual encoding
---------------
- Modern reference (718): uniform light-grey semi-transparent background.
- Ancient samples:
    * SHAPE = core (field_sample_id): CAM2509=circle, CAM23-13=square,
      CAM23-11=diamond, CAM22-08=triangle-up, CAM2201=x.
    * COLOR = continuous Age CE (default 970-2021) via a Viridis colorscale.
    * Hover: sample_id / core / site / depth_cm / age_CE / libtype / prep / besthit.

Views (dropdown)
----------------
1. Per-core isolation : show one core, hide others.
2. Panoramic + age slider : all cores, range slider on Age CE filters points.
3. Trajectory : per-core centroid path over 100-yr windows (PC1-PC2 plane,
   no-arrow polyline by default, optional arrow toggle) + PC-vs-Age curve.
4. QC sub-library comparison : connect each base_robot's multi-lib points
   (SG / C1 / C2 / BH / pooled) with dashed lines.

Static exports (if matplotlib available)
----------------------------------------
- {prefix}_panoramic.png, {prefix}_percore.png
"""
import argparse
import collections
import sys
from pathlib import Path


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--evec", required=True)
    ap.add_argument("--eval", required=True)
    ap.add_argument("--meta", required=True,
                    help="enriched 440 meta TSV: sample_id base_robot core site "
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


def collect(args):
    evals = load_eval(args.eval)
    total = sum(evals) or 1.0
    meta = load_meta(args.meta)
    rows = load_evec(args.evec)
    modern = [r for r in rows if r[2] != "Ancient"]
    anc_rows = [r for r in rows if r[2] == "Ancient"]

    ancient = collections.defaultdict(list)
    by_id = {}
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
        by_id[iid] = rec

    cores = [c for c in CORE_ORDER if c in ancient]
    return evals, total, modern, ancient, by_id, cores


def _trace_visible(base, per_core):
    """visible list for the 2 main views: index0=modern, then per-core traces."""
    return [base] + per_core


def trajectory(ancient, core, age_min, age_max, bin=100):
    """Return [(age_label, x_centroid, y_centroid)] for one core over 100-yr bins."""
    pts = [p for p in ancient[core] if p["age"] is not None]
    bins = collections.defaultdict(list)
    for p in pts:
        b = int(p["age"] // bin) * bin
        bins[b].append(p)
    out = []
    for b in sorted(bins):
        xs = [p["x"] for p in bins[b]]
        ys = [p["y"] for p in bins[b]]
        out.append((f"{b}s", sum(xs) / len(xs), sum(ys) / len(ys)))
    return out


def build_html(args, evals, total, modern, ancient, cores):
    import plotly.graph_objects as go

    # trace 0: modern background
    traces = [go.Scatter(
        x=[p[1][0] for p in modern], y=[p[1][1] for p in modern],
        mode="markers", name="Modern (718)",
        marker=dict(color="#C0C0C0", size=6, opacity=0.35), hoverinfo="skip",
        showlegend=False)]
    core_trace_idx = {}
    for core in cores:
        pts = ancient[core]
        if not pts:
            continue
        tr = go.Scatter(
            x=[p["x"] for p in pts], y=[p["y"] for p in pts], mode="markers",
            name=core,
            marker=dict(size=9, symbol=CORE_SYMBOL.get(core, "circle"),
                        color=[p["age"] for p in pts],
                        colorscale="Viridis", cmin=args.age_min, cmax=args.age_max,
                        colorbar=dict(title="Age CE"),
                        line=dict(width=0.5, color="black")),
            customdata=[p["hover"] for p in pts],
            hovertemplate="%{customdata}<extra></extra>", showlegend=True)
        core_trace_idx[core] = len(traces)
        traces.append(tr)
    n_core = len(core_trace_idx)

    # ---- view 3: trajectory traces (polyline through 100-yr centroids) ----
    traj_traces = []
    traj_core_idx = {}
    for core in cores:
        tr = trajectory(ancient, core, args.age_min, args.age_max)
        if len(tr) < 2:
            continue
        xs = [p[1] for p in tr]
        ys = [p[2] for p in tr]
        labs = [p[0] for p in tr]
        t = go.Scatter(
            x=xs, y=ys, mode="lines+markers+text", name=f"{core} traj",
            line=dict(color="black", width=2), marker=dict(size=9, color="red"),
            text=labs, textposition="bottom right", textfont=dict(size=9),
            hoverinfo="skip", showlegend=False)
        traj_core_idx[core] = len(traces) + len(traj_traces)
        traj_traces.append(t)

    # ---- view 4: QC sub-library comparison traces ----
    libgroup = collections.defaultdict(list)
    for core in cores:
        for p in ancient[core]:
            if p["libtype"] in ("SG", "C1", "C2", "pooled_nobesthit", "pooled_besthit"):
                libgroup[p["base"]].append(p)
    qc_traces = []
    order = {"SG": 0, "C1": 1, "C2": 2, "pooled_nobesthit": 3, "pooled_besthit": 4}
    for base, pts in libgroup.items():
        if len(pts) < 2:
            continue
        pts = sorted(pts, key=lambda p: order.get(p["libtype"], 9))
        xs = [p["x"] for p in pts] + [pts[0]["x"]]
        ys = [p["y"] for p in pts] + [pts[0]["y"]]
        txt = [p["libtype"] for p in pts] + [pts[0]["libtype"]]
        qc_traces.append(go.Scatter(
            x=xs, y=ys, mode="lines+markers+text", name=f"{base} QC",
            line=dict(color="rgba(0,0,0,0.6)", width=1, dash="dot"),
            marker=dict(size=7, color="black"),
            text=txt, textposition="top right", textfont=dict(size=8),
            customdata=[p["hover"] for p in pts] + [pts[0]["hover"]],
            hovertemplate="%{customdata}<extra></extra>", showlegend=False))

    n_traj = len(traj_traces)
    n_qc = len(qc_traces)
    all_traces = traces + traj_traces + qc_traces

    # visibility helpers (over all_traces)
    core_vis = lambda mask: [True] + mask + [0] * n_traj + [0] * n_qc
    qc_off = [0] * n_qc

    # ---- view 1: per-core isolation ----
    buttons1 = [dict(label="All cores", method="update",
                     args=[{"visible": [True] + [1] * n_core + [0] * n_traj + qc_off}])]
    for core in cores:
        if core not in core_trace_idx:
            continue
        vis = [True] + [1 if c == core else 0 for c in cores] + [0] * n_traj + qc_off
        buttons1.append(dict(label=f"Core: {core}", method="update",
                             args=[{"visible": vis}]))

    # ---- view 2: age slider (bin buttons) ----
    age_buttons = [dict(label="All ages", method="update",
                        args=[{"visible": [True] + [1] * n_core + [0] * n_traj + qc_off}])]
    for lo in range(int(args.age_min), int(args.age_max), 100):
        hi = min(lo + 100, int(args.age_max))
        vis = [True]
        for core in cores:
            a = [p["age"] for p in ancient[core]]
            vis.append([1 if (v is not None and lo <= v < hi) else 0 for v in a])
        vis += [0] * n_traj + qc_off
        age_buttons.append(dict(label=f"{lo}-{hi} CE", method="update",
                                args=[{"visible": vis}]))

    # ---- view 3: trajectory button ----
    vis3 = [True] + [1] * n_core + [1] * n_traj + qc_off

    # ---- view 4: QC button ----
    vis4 = [True] + [1] * n_core + [0] * n_traj + [1] * n_qc

    view_buttons = [
        dict(label="Panoramic (all ages)", method="update",
             args=[{"visible": [True] + [1] * n_core + [0] * n_traj + qc_off}]),
        dict(label="Trajectory (100-yr)", method="update", args=[{"visible": vis3}]),
        dict(label="QC sub-library", method="update", args=[{"visible": vis4}]),
    ]

    fig = go.Figure(data=all_traces)
    fig.update_layout(
        title=dict(text=f"{args.title} | markers={args.nmarkers}", font=dict(size=16)),
        hovermode="closest",
        xaxis=dict(title=f"PC1 ({evals[0]/total*100:.2f}%)"),
        yaxis=dict(title=f"PC2 ({evals[1]/total*100:.2f}%)"),
        height=780, width=1050,
        margin=dict(l=60, r=60, t=140, b=100),
        updatemenus=[
            dict(buttons=view_buttons, direction="down", showactive=True,
                 x=0.0, y=1.24, xanchor="left", yanchor="top",
                 bgcolor="rgba(240,240,240,0.9)"),
            dict(buttons=buttons1, direction="down", showactive=True,
                 x=0.0, y=1.02, xanchor="left", yanchor="top",
                 bgcolor="rgba(240,240,240,0.9)"),
            dict(buttons=age_buttons, direction="up", showactive=True,
                 x=1.0, y=-0.15, xanchor="right", yanchor="top",
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

    def draw(ax, data):
        ax.scatter([p[1][0] for p in modern], [p[1][1] for p in modern],
                   s=8, color="#C0C0C0", alpha=0.35, linewidths=0)
        for c, pts in data.items():
            ages = [p["age"] for p in pts]
            ax.scatter([p["x"] for p in pts], [p["y"] for p in pts], c=ages,
                       cmap="viridis", s=40, vmin=args.age_min, vmax=args.age_max,
                       edgecolor="black", linewidth=0.5)
        ax.set_xlabel(f"PC1 ({evals[0]/total*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({evals[1]/total*100:.1f}%)")

    # panoramic
    fig, ax = plt.subplots(figsize=(8, 7))
    draw(ax, dict(ancient))
    sm = plt.cm.ScalarMappable(cmap="viridis",
                               norm=plt.Normalize(args.age_min, args.age_max))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Age CE")
    ax.set_title(f"{args.title} | markers={args.nmarkers}")
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_panoramic.png", dpi=150)
    print(f"wrote {args.out_prefix}_panoramic.png")

    # per-core facets
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

    # PC-vs-Age curve (PC1 per core over Age CE)
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
    evals, total, modern, ancient, by_id, cores = collect(args)
    build_html(args, evals, total, modern, ancient, cores)
    if not args.no_png:
        try:
            build_pngs(args, evals, total, modern, ancient, cores)
        except ImportError as e:
            print(f"WARNING: matplotlib not available, skipping PNG ({e})",
                  file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
