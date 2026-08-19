#!/usr/bin/env python3
"""Plot a smartpca .evec/.eval pair.

Modern samples are colored by population label with semi-transparent markers so
overlapping points show through. Ancient samples are hollow black circles; a
configurable set of low-coverage/unreliable ancient samples is drawn as hollow
triangles. Labels default to the low-confidence samples only (no 16-label pileup);
pass --label-all to annotate every ancient sample, using adjustText for de-overlap
when the library is available. One shared legend is placed outside the figure.
Axis labels show per-PC explained variance; the title shows the marker count.
"""

import argparse
import sys
from collections import defaultdict


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evec", required=True)
    ap.add_argument("--eval", required=True)
    ap.add_argument("--ind", required=True,
                    help="accepted for CLI compatibility; labels are read from "
                         "the .evec's last column instead")
    ap.add_argument("--nmarkers", type=int, required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--lowconf-samples", default="LV7008416294",
                    help="comma-separated ancient IDs drawn as hollow triangles")
    ap.add_argument("--label-all", action="store_true",
                    help="label every ancient sample (uses adjustText if installed)")
    ap.add_argument("--modern-alpha", type=float, default=0.5)
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


def main():
    args = parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lowconf = {s.strip() for s in args.lowconf_samples.split(",") if s.strip()}

    evals = load_eval(args.eval)
    total = sum(evals) or 1.0
    rows = load_evec(args.evec)
    npc = len(rows[0][1])
    if len(evals) <= npc:
        print(f"WARNING: {args.eval} has only {len(evals)} eigenvalues (<= {npc} PCs). "
              f"PC% is normalized over the top-{len(evals)} only and will be inflated. "
              f"smartpca .eval should contain the FULL spectrum (sum = n_samples - 1).",
              file=sys.stderr)

    by_lab = defaultdict(list)
    ancient = []
    for iid, vals, label in rows:
        if label == "Ancient":
            ancient.append((iid, vals))
        else:
            by_lab[label].append((iid, vals))

    normal_ancient = [(i, v) for i, v in ancient if i not in lowconf]
    lowconf_ancient = [(i, v) for i, v in ancient if i in lowconf]

    npairs = min(5, npc // 2)
    fig, axes = plt.subplots(1, npairs, figsize=(6.8 * npairs, 5.9), squeeze=False)
    cmap = plt.get_cmap("tab20")
    labs = sorted(by_lab)

    # Build one shared legend (placed outside the axes) instead of per-panel legends.
    handles, labels = [], []
    for li, lab in enumerate(labs):
        handles.append(plt.Line2D([0], [0], marker="o", color="none",
                                  markerfacecolor=cmap(li % 20), markeredgecolor="none",
                                  markersize=6, alpha=args.modern_alpha))
        labels.append(lab)
    if normal_ancient:
        handles.append(plt.Line2D([0], [0], marker="o", color="none",
                                  markerfacecolor="none", markeredgecolor="black",
                                  markersize=7, linewidth=1.2))
        labels.append(f"Ancient (n={len(normal_ancient)})")
    if lowconf_ancient:
        handles.append(plt.Line2D([0], [0], marker="^", color="none",
                                  markerfacecolor="none", markeredgecolor="black",
                                  markersize=8, linewidth=1.4))
        labels.append("Ancient low-conf (" + ",".join(sorted(lowconf)) + ")")

    to_label = ancient if args.label_all else lowconf_ancient
    panel_texts = [[] for _ in range(npairs)]

    for pi in range(npairs):
        ax = axes[0][pi]
        xi, yi = 2 * pi, 2 * pi + 1
        for li, lab in enumerate(labs):
            xs = [v[xi] for _, v in by_lab[lab]]
            ys = [v[yi] for _, v in by_lab[lab]]
            ax.scatter(xs, ys, s=args.modern_size, color=cmap(li % 20),
                       alpha=args.modern_alpha, linewidths=0, zorder=1)
        if normal_ancient:
            ax.scatter([v[xi] for _, v in normal_ancient], [v[yi] for _, v in normal_ancient],
                       s=30, marker="o", facecolor="none", edgecolor="black",
                       linewidth=1.2, zorder=3)
        if lowconf_ancient:
            ax.scatter([v[xi] for _, v in lowconf_ancient], [v[yi] for _, v in lowconf_ancient],
                       s=42, marker="^", facecolor="none", edgecolor="black",
                       linewidth=1.4, zorder=4)
        for iid, v in to_label:
            t = ax.annotate(iid, (v[xi], v[yi]), textcoords="offset points",
                            xytext=(4, 4), fontsize=7, fontweight="bold", zorder=5)
            panel_texts[pi].append(t)
        ax.set_xlabel(f"PC{xi + 1} ({evals[xi] / total * 100:.2f}%)")
        ax.set_ylabel(f"PC{yi + 1} ({evals[yi] / total * 100:.2f}%)")
        ax.set_title(f"{args.title} | markers={args.nmarkers}")

    if args.label_all:
        try:
            from adjustText import adjust_text
        except ImportError:
            adjust_text = None
        if adjust_text is not None:
            for pi in range(npairs):
                if panel_texts[pi]:
                    adjust_text(panel_texts[pi], ax=axes[0][pi], expand=(1.2, 1.4))

    fig.legend(handles=handles, labels=labels, loc="center left",
               bbox_to_anchor=(1.02, 0.5), fontsize=7, markerscale=1.2)
    fig.tight_layout(rect=[0, 0, 0.86, 1])
    out = f"{args.out_prefix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
