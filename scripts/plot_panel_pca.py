#!/usr/bin/env python3
"""Plot modern-only PCA from plink2 --pca output, colored by .ind population
label, with per-axis explained variance and marker count in the title."""
import argparse
from collections import defaultdict


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eigenvec", required=True)
    ap.add_argument("--eigenval", required=True)
    ap.add_argument("--ind", required=True)
    ap.add_argument("--nmarkers", type=int, required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--out-prefix", required=True)
    return ap.parse_args()


def load_ind(path):
    lab = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if not f:
                continue
            lab[f[0]] = f[2] if len(f) >= 3 else "?"
    return lab


def load_eigenvec(path):
    rows = []
    with open(path) as fh:
        head = fh.readline().split()
        npc = len(head) - 2
        for line in fh:
            f = line.split()
            if len(f) < 2 + npc:
                continue
            rows.append((f[1], [float(x) for x in f[2:2 + npc]]))
    return rows, npc


def load_eigenval(path):
    vals = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                vals.append(float(line))
    return vals


def main():
    args = parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = load_ind(args.ind)
    rows, npc = load_eigenvec(args.eigenvec)
    evals = load_eigenval(args.eigenval)
    total = sum(evals) or 1.0
    if len(evals) <= npc:
        import sys as _sys
        print(f"WARNING: {args.eigenval} has only {len(evals)} eigenvalues (<= {npc} PCs). "
              f"PC% is normalized over the top-{len(evals)} only and will be inflated. "
              f"Re-run `plink2 --pca` WITHOUT a count so the full spectrum is written.",
              file=_sys.stderr)
    by_lab = defaultdict(list)
    for iid, vals in rows:
        by_lab[labels.get(iid, "?")].append(vals)

    npairs = min(5, npc // 2)
    fig, axes = plt.subplots(1, npairs, figsize=(6.8 * npairs, 5.8), squeeze=False)
    cmap = plt.get_cmap("tab20")
    labs = sorted(by_lab)
    for pi in range(npairs):
        ax = axes[0][pi]
        xi, yi = 2 * pi, 2 * pi + 1
        for li, lab in enumerate(labs):
            xs = [v[xi] for v in by_lab[lab]]
            ys = [v[yi] for v in by_lab[lab]]
            ax.scatter(xs, ys, s=12, color=cmap(li % 20), label=lab, alpha=0.7, linewidths=0)
        ax.set_xlabel(f"PC{xi + 1} ({evals[xi] / total * 100:.2f}%)")
        ax.set_ylabel(f"PC{yi + 1} ({evals[yi] / total * 100:.2f}%)")
        ax.set_title(f"{args.title} | markers={args.nmarkers}")
        ax.legend(fontsize=6, markerscale=1.5, loc="best")
    fig.tight_layout()
    out = f"{args.out_prefix}.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
