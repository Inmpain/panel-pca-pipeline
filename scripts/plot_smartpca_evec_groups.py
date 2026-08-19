#!/usr/bin/env python3
"""Plot a smartpca .evec/.eval pair, coloring ancient samples by group (e.g.
angkor vs nanzuo) with a shared legend, suptitle, and a 2x3 subplot grid.

Modern samples are semi-transparent and colored by population label. Ancient
samples are drawn as larger solid circles (or hollow triangles for
--lowconf-samples) colored by their group, standing out against the modern
background. Text labels are OFF by default (they clutter); enable with
--label-ancient (all) or --label-lowconf (only low-confidence samples).

--ancient-meta is a 3-column TSV: sample_id<TAB>group<TAB>label
  e.g.
    LV7008416379<TAB>angkor<TAB>53_1709
    YWL1-A3483<TAB>nanzuo<TAB>YWL1-A3483
"""
import argparse
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
                    help="TSV: sample_id<TAB>group<TAB>label (one ancient per line)")
    ap.add_argument("--lowconf-samples", default="",
                    help="comma-separated ancient IDs drawn as hollow triangles")
    ap.add_argument("--label-ancient", action="store_true",
                    help="label every ancient sample (default off)")
    ap.add_argument("--label-lowconf", action="store_true",
                    help="label only low-confidence ancient samples (default off)")
    ap.add_argument("--modern-alpha", type=float, default=0.45)
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
            meta[sid] = (group, label)
    return meta


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
    ancient = []  # (iid, vals, group, label, is_lowconf)
    for iid, vals, label in rows:
        if label == "Ancient":
            group, alabel = meta.get(iid, ("?", iid))
            ancient.append((iid, vals, group, alabel, iid in lowconf))
        else:
            by_lab[label].append((iid, vals))

    groups = sorted({g for _, _, g, _, _ in ancient})
    group_cmap = {g: plt.get_cmap("tab10")(i % 10) for i, g in enumerate(groups)}

    npairs = min(5, npc // 2)
    ncols, nrows = 3, 2  # 5 pairs in a 2x3 grid
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.0 * ncols, 5.0 * nrows),
                             squeeze=False)
    cmap = plt.get_cmap("tab20")
    labs = sorted(by_lab)

    # --- legend handles: modern populations + ancient groups ---
    handles, labels = [], []
    for li, lab in enumerate(labs):
        handles.append(plt.Line2D([0], [0], marker="o", color="none",
                                  markerfacecolor=cmap(li % 20), markeredgecolor="none",
                                  markersize=8, alpha=args.modern_alpha))
        labels.append(lab)
    for g in groups:
        handles.append(plt.Line2D([0], [0], marker="o", color="none",
                                  markerfacecolor=group_cmap[g], markeredgecolor="black",
                                  markersize=10, markeredgewidth=1.3))
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
        for iid, v, g, alabel, is_low in ancient:
            if is_low:
                ax.scatter(v[xi], v[yi], s=110, marker="^", facecolor="none",
                           edgecolor=group_cmap[g], linewidth=2.2, zorder=4)
            else:
                ax.scatter(v[xi], v[yi], s=80, marker="o", color=group_cmap[g],
                           edgecolor="black", linewidth=1.3, zorder=3)
            do_label = args.label_ancient or (is_low and args.label_lowconf)
            if do_label:
                t = ax.annotate(alabel, (v[xi], v[yi]), textcoords="offset points",
                                xytext=(5, 5), fontsize=7, fontweight="bold", zorder=5)
                panel_texts[pi].append(t)
        ax.set_xlabel(f"PC{xi + 1} ({evals[xi] / total * 100:.2f}%)", fontsize=10)
        ax.set_ylabel(f"PC{yi + 1} ({evals[yi] / total * 100:.2f}%)", fontsize=10)
        ax.set_title(f"PC{xi + 1} vs PC{yi + 1}", fontsize=11)

    # hide any unused cells (the 6th in a 2x3 grid)
    for pi in range(npairs, nrows * ncols):
        axes[pi // ncols][pi % ncols].axis("off")

    # de-overlap labels only if any are drawn
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

    fig.suptitle(f"{args.title} | markers={args.nmarkers}", fontsize=14, y=0.995)
    fig.legend(handles=handles, labels=labels, loc="center left",
               bbox_to_anchor=(1.02, 0.5), fontsize=9, markerscale=1.4,
               frameon=False)
    fig.tight_layout(rect=[0, 0, 0.85, 0.95])
    out = f"{args.out_prefix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
