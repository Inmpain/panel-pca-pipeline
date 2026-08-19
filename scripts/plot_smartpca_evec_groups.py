#!/usr/bin/env python3
"""Plot a smartpca .evec/.eval pair, coloring ancient samples by group (e.g.
angkor vs nanzuo), output both a static PNG (matplotlib 2x3 overview) and an
interactive single-panel HTML (plotly).

HTML: one large panel at a time, dropdown switches PC pairs. A second dropdown
toggles "arrows + depth labels": clicking it reveals per-group depth-ordered
trajectory arrows (shortened so they don't pierce points) plus a depth text
label next to each angkor triangle (e.g. "53cm"). Modern samples are low-opacity
background; ancient samples are high-saturation larger markers.

--ancient-meta is a TSV: sample_id<TAB>group<TAB>label[<TAB>order]
  order = numeric ordering for arrows (descending = deep/old -> shallow/new).
  Depth text is auto-derived when the label starts with "depth_age" (e.g. "53_1709"
  -> "53cm"); otherwise no label.
"""
import argparse
import sys
from collections import defaultdict


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evec", required=True)
    ap.add_argument("--eval", required=True)
    ap.add_argument("--ind", required=True,
                    help="accepted for CLI compatibility; labels are read from .evec last column")
    ap.add_argument("--nmarkers", type=int, required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--ancient-meta", required=True,
                    help="TSV: sample_id<TAB>group<TAB>label[<TAB>order]")
    ap.add_argument("--lowconf-samples", default="",
                    help="comma-separated ancient IDs drawn as hollow triangles")
    ap.add_argument("--label-ancient", action="store_true",
                    help="label every ancient sample in the PNG (default off)")
    ap.add_argument("--label-lowconf", action="store_true",
                    help="label only low-confidence ancient samples in the PNG (default off)")
    ap.add_argument("--no-html", action="store_true",
                    help="skip the plotly HTML export")
    ap.add_argument("--modern-alpha", type=float, default=0.30)
    ap.add_argument("--modern-size", type=float, default=12.0)
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
            iid = f[0]
            label = f[-1]
            vals = list(map(float, f[1:-1]))
            rows.append((iid, vals, label))
    return rows


def load_meta(path):
    """Return {sample_id: (group, label, order_or_None)}."""
    meta = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 3:
                continue
            sid, group, label = f[0], f[1], f[2]
            order = None
            if len(f) >= 4:
                try:
                    order = float(f[3])
                except ValueError:
                    order = None
            meta[sid] = (group, label, order)
    return meta


def depth_text(label):
    """Auto depth label from 'depth_age' style labels ('53_1709' -> '53cm')."""
    if label and label[0].isdigit() and "_" in label:
        return label.split("_", 1)[0] + "cm"
    return ""


def pair_coords(by_lab, ancient, xi, yi):
    xs, ys = [], []
    for lab in by_lab:
        for _, v in by_lab[lab]:
            xs.append(v[xi])
            ys.append(v[yi])
    for _, v, _, _, _, _ in ancient:
        xs.append(v[xi])
        ys.append(v[yi])
    return xs, ys


def pct_limits(vals):
    import numpy as np
    a = np.asarray(vals, dtype=float)
    lo, hi = np.nanpercentile(a, 1.0), np.nanpercentile(a, 99.0)
    pad = (hi - lo) * 0.12
    if not np.isfinite(pad) or pad <= 0:
        pad = 1.0
    return float(lo - pad), float(hi + pad)


def export_html(args, evals, total, by_lab, ancient, groups, group_cmap, labs,
                npairs):
    """Single-panel HTML: PC-pair dropdown + optional arrows/depth labels."""
    try:
        import matplotlib
        import matplotlib.colors as mcolors
        import plotly.graph_objects as go
    except ImportError:
        print("WARNING: plotly (or matplotlib) not installed, skipping HTML export",
              file=sys.stderr)
        return

    cmap = matplotlib.cm.get_cmap("tab20")
    lab_color = {lab: mcolors.to_hex(cmap(i % 20)) for i, lab in enumerate(labs)}
    group_color = {g: mcolors.to_hex(group_cmap[g]) for g in groups}

    group_pts = {g: [] for g in groups}
    for iid, v, g, alabel, is_low, order in ancient:
        if g in group_pts:
            group_pts[g].append((order, iid, v, alabel, is_low))

    traces = []
    trace_pairs = []
    trace_is_ancient = []
    annotations = []
    ann_pairs = []

    for pi in range(npairs):
        xi, yi = 2 * pi, 2 * pi + 1
        ht = (f"%{{customdata}}<br>PC{xi + 1}=%{{x:.4f}}<br>"
              f"PC{yi + 1}=%{{y:.4f}}<extra></extra>")

        for lab in labs:
            xs = [v[xi] for _, v in by_lab[lab]]
            ys = [v[yi] for _, v in by_lab[lab]]
            cd = [f"{iid} ({lab})" for iid, _ in by_lab[lab]]
            traces.append(go.Scatter(
                x=xs, y=ys, mode="markers", name=lab,
                marker=dict(color=lab_color[lab], size=6, opacity=0.4),
                customdata=cd, hovertemplate=ht, showlegend=True,
            ))
            trace_pairs.append(pi)
            trace_is_ancient.append(False)

        for g in groups:
            pts = group_pts[g]
            xs = [p[2][xi] for p in pts]
            ys = [p[2][yi] for p in pts]
            cd = [f"{p[1]} · {p[3]} · {g}" for p in pts]
            txt = [depth_text(p[3]) for p in pts]
            sym = "triangle-up" if any(p[4] for p in pts) else "circle"
            traces.append(go.Scatter(
                x=xs, y=ys, mode="markers", name=g,
                marker=dict(color=group_color[g], size=8, symbol=sym, opacity=0.9,
                            line=dict(color="black", width=0.8)),
                customdata=cd, hovertemplate=ht, showlegend=True,
                text=txt, textposition="top right",
                textfont=dict(size=9, color="black"),
            ))
            trace_pairs.append(pi)
            trace_is_ancient.append(True)

        # arrows: descending order (deep/old -> shallow/new), shortened 8%
        for g in groups:
            ordered = sorted([p for p in group_pts[g] if p[0] is not None],
                             key=lambda t: t[0], reverse=True)
            for k in range(len(ordered) - 1):
                tail = ordered[k][2]
                head = ordered[k + 1][2]
                dx = head[xi] - tail[xi]
                dy = head[yi] - tail[yi]
                sx = tail[xi] + 0.08 * dx
                sy = tail[yi] + 0.08 * dy
                ex = head[xi] - 0.08 * dx
                ey = head[yi] - 0.08 * dy
                annotations.append(dict(
                    x=ex, y=ey, ax=sx, ay=sy,
                    xref="x", yref="y", axref="x", ayref="y",
                    showarrow=True, arrowhead=2, arrowsize=1.1, arrowwidth=2.0,
                    arrowcolor=group_color[g], opacity=0.8, visible=True,
                ))
                ann_pairs.append(pi)

    fig = go.Figure()
    for t in traces:
        fig.add_trace(t)

    n = len(traces)
    init_vis = [p == 0 for p in trace_pairs]
    for k, t in enumerate(fig.data):
        t.visible = init_vis[k]
    fig.update_layout(annotations=[])  # no arrows by default

    def mode_array(with_text):
        out = []
        for k, t in enumerate(traces):
            if trace_is_ancient[k]:
                out.append("markers+text" if with_text else "markers")
            else:
                out.append(t.mode)
        return out

    buttons_plain = []
    buttons_arrow = []
    for pi in range(npairs):
        xi, yi = 2 * pi, 2 * pi + 1
        vis = [p == pi for p in trace_pairs]
        anns = [dict(a) for a, ap in zip(annotations, ann_pairs) if ap == pi]
        xs_all, ys_all = pair_coords(by_lab, ancient, xi, yi)
        layout = {
            "xaxis.title.text": f"PC{xi + 1} ({evals[xi] / total * 100:.2f}%)",
            "yaxis.title.text": f"PC{yi + 1} ({evals[yi] / total * 100:.2f}%)",
            "xaxis.range": pct_limits(xs_all),
            "yaxis.range": pct_limits(ys_all),
        }
        lbl = f"PC{xi + 1} vs PC{yi + 1}"
        buttons_plain.append(dict(
            label=lbl, method="update",
            args=[{"visible": vis, "mode": mode_array(False)},
                  {**layout, "annotations": [dict(a, visible=False) for a in anns]}]))
        buttons_arrow.append(dict(
            label=lbl, method="update",
            args=[{"visible": vis, "mode": mode_array(True)},
                  {**layout, "annotations": [dict(a, visible=True) for a in anns]}]))

    fig.update_layout(
        title=dict(text=f"{args.title} | markers={args.nmarkers}", font=dict(size=16)),
        height=820, width=1150,
        updatemenus=[
            dict(buttons=buttons_plain, direction="down", showactive=True,
                 x=0.02, y=1.16, xanchor="left", yanchor="top",
                 font=dict(size=12), bgcolor="rgba(240,240,240,0.9)",
                 bordercolor="#ccc", borderwidth=1,
                 pad=dict(r=6)),
            dict(buttons=buttons_arrow, direction="down", showactive=True,
                 x=0.16, y=1.16, xanchor="left", yanchor="top",
                 font=dict(size=12), bgcolor="rgba(210,240,210,0.9)",
                 bordercolor="#8c8", borderwidth=1,
                 pad=dict(r=6)),
        ],
        legend=dict(font=dict(size=12), x=1.02, y=1.0, xanchor="left", yanchor="top"),
        hovermode="closest",
        margin=dict(l=70, r=220, t=100, b=60),
    )
    html_out = f"{args.out_prefix}.html"
    fig.write_html(html_out)
    print(f"wrote {html_out}")


def main():
    args = parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    meta = load_meta(args.ancient_meta)
    lowconf = {s.strip() for s in args.lowconf_samples.split(",") if s.strip()}

    evals = load_eval(args.eval)
    total = sum(evals) or 1.0
    rows = load_evec(args.evec)
    npc = len(rows[0][1])

    by_lab = defaultdict(list)
    ancient = []  # (iid, vals, group, label, is_lowconf, order)
    for iid, vals, label in rows:
        if label == "Ancient":
            group, alabel, order = meta.get(iid, ("?", iid, None))
            ancient.append((iid, vals, group, alabel, iid in lowconf, order))
        else:
            by_lab[label].append((iid, vals))

    groups = sorted({g for _, _, g, _, _, _ in ancient})
    group_cmap = {g: plt.get_cmap("tab10")(i % 10) for i, g in enumerate(groups)}

    npairs = min(5, npc // 2)
    ncols, nrows = 3, 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.0 * ncols, 5.0 * nrows),
                             squeeze=False)
    cmap = plt.get_cmap("tab20")
    labs = sorted(by_lab)

    handles, labels = [], []
    for li, lab in enumerate(labs):
        handles.append(plt.Line2D([0], [0], marker="o", color="none",
                                  markerfacecolor=cmap(li % 20), markeredgecolor="none",
                                  markersize=8, alpha=args.modern_alpha))
        labels.append(lab)
    for g in groups:
        handles.append(plt.Line2D([0], [0], marker="o", color="none",
                                  markerfacecolor=group_cmap[g], markeredgecolor="black",
                                  markersize=10, markeredgewidth=1.0))
        labels.append(g)

    panel_texts = [[] for _ in range(npairs)]
    for pi in range(npairs):
        ax = axes[pi // ncols][pi % ncols]
        xi, yi = 2 * pi, 2 * pi + 1
        for li, lab in enumerate(labs):
            xs = [v[xi] for _, v in by_lab[lab]]
            ys = [v[yi] for _, v in by_lab[lab]]
            ax.scatter(xs, ys, s=args.modern_size, color=cmap(li % 20),
                       alpha=args.modern_alpha, linewidths=0, zorder=1)
        for iid, v, g, alabel, is_low, _order in ancient:
            if is_low:
                ax.scatter(v[xi], v[yi], s=70, marker="^", facecolor="none",
                           edgecolor=group_cmap[g], linewidth=1.4, zorder=4)
            else:
                ax.scatter(v[xi], v[yi], s=60, marker="o", color=group_cmap[g],
                           edgecolor="black", linewidth=1.0, zorder=3)
            do_label = args.label_ancient or (is_low and args.label_lowconf)
            if do_label:
                t = ax.annotate(alabel, (v[xi], v[yi]), textcoords="offset points",
                                xytext=(5, 5), fontsize=7, fontweight="bold", zorder=5)
                panel_texts[pi].append(t)

        xs_all, ys_all = pair_coords(by_lab, ancient, xi, yi)
        ax.set_xlim(*pct_limits(xs_all))
        ax.set_ylim(*pct_limits(ys_all))

        ax.set_xlabel(f"PC{xi + 1} ({evals[xi] / total * 100:.2f}%)", fontsize=10)
        ax.set_ylabel(f"PC{yi + 1} ({evals[yi] / total * 100:.2f}%)", fontsize=10)
        ax.set_title(f"PC{xi + 1} vs PC{yi + 1}", fontsize=11)

    for pi in range(npairs, nrows * ncols):
        axes[pi // ncols][pi % ncols].axis("off")

    if args.label_ancient or args.label_lowconf:
        try:
            from adjustText import adjust_text
        except ImportError:
            adjust_text = None
        if adjust_text is not None:
            for pi in range(npairs):
                if panel_texts[pi]:
                    adjust_text(panel_texts[pi], ax=axes[pi // ncols][pi % ncols],
                                expand=(1.2, 1.4))

    fig.suptitle(f"{args.title} | markers={args.nmarkers}", fontsize=15,
                 y=0.98, fontweight="bold")
    fig.legend(handles=handles, labels=labels, loc="center left",
               bbox_to_anchor=(1.02, 0.5), fontsize=9, markerscale=1.4,
               frameon=False)
    fig.tight_layout(rect=[0, 0, 0.85, 0.94])
    out = f"{args.out_prefix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", pad_inches=0.3)
    print(f"wrote {out}")

    if not args.no_html:
        export_html(args, evals, total, by_lab, ancient, groups, group_cmap,
                    labs, npairs)


if __name__ == "__main__":
    main()
