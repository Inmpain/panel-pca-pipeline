#!/usr/bin/env python3
"""Ancient-sample -> modern-population centroid distance ranking in PCA space,
with bootstrap stability of the nearest assignment. stdlib only.

Usage: python3 distance_ranking.py EVEC OUT_PREFIX [NPC] [N_BOOT]
  NPC    = number of PCs to use (default 2 = PC1-PC2)
  N_BOOT = bootstrap resamples (default 200)

Outputs:
  OUT_PREFIX.summary.tsv : nearest / 2nd nearest / gap / d1_over_d2 / n_nearest /
                           n_pops / nearest_boot_pct
  OUT_PREFIX.top5.tsv    : top-5 populations per sample with distances
"""
import sys
import random
from collections import defaultdict


def parse(evec, npc):
    modern = defaultdict(list)
    ancient = {}
    for line in open(evec):
        if line.lstrip().startswith("#"):
            continue
        f = line.split()
        if len(f) < npc + 2:
            continue
        sid, lab = f[0], f[-1]
        c = [float(x) for x in f[1:npc + 1]]
        if lab == "Ancient":
            ancient[sid] = c
        else:
            modern[lab].append(c)
    return modern, ancient


def centroid(pts, npc):
    return [sum(p[i] for p in pts) / len(pts) for i in range(npc)]


def dist(a, b):
    return sum((a[i] - b[i]) ** 2 for i in range(len(a))) ** 0.5


def main():
    evec, prefix = sys.argv[1], sys.argv[2]
    npc = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    nboot = int(sys.argv[4]) if len(sys.argv) > 4 else 200
    random.seed(42)

    modern, ancient = parse(evec, npc)
    pops = sorted(modern)
    print("populations:", {p: len(modern[p]) for p in pops})
    if not ancient:
        print("FATAL: no Ancient-labeled rows in evec")
        sys.exit(1)

    cent = {p: centroid(modern[p], npc) for p in pops}

    boot = {s: defaultdict(int) for s in ancient}
    for _ in range(nboot):
        bc = {p: centroid([random.choice(modern[p]) for _ in modern[p]], npc) for p in pops}
        for s, c in ancient.items():
            best = min(pops, key=lambda p: dist(c, bc[p]))
            boot[s][best] += 1

    with open(prefix + ".summary.tsv", "w") as o:
        o.write("sample\tnearest_pop\tnearest_d\t2nd_pop\t2nd_d\tgap_d2d1\td1_over_d2\t"
                "n_nearest\tn_pops\tnearest_boot_pct\n")
        for s in sorted(ancient):
            c = ancient[s]
            r = sorted(pops, key=lambda p: dist(c, cent[p]))
            d1, d2 = dist(c, cent[r[0]]), dist(c, cent[r[1]])
            frac = boot[s][r[0]] / nboot
            o.write(f"{s}\t{r[0]}\t{d1:.5f}\t{r[1]}\t{d2:.5f}\t{d2 - d1:.5f}\t"
                    f"{d1 / d2:.3f}\t{len(modern[r[0]])}\t{len(pops)}\t{frac:.3f}\n")

    with open(prefix + ".top5.tsv", "w") as o:
        o.write("sample\trank\tpop\tdistance\tn_pop\n")
        for s in sorted(ancient):
            c = ancient[s]
            r = sorted(pops, key=lambda p: dist(c, cent[p]))
            for rank, p in enumerate(r[:5], 1):
                o.write(f"{s}\t{rank}\t{p}\t{dist(c, cent[p]):.5f}\t{len(modern[p])}\n")

    print("wrote", prefix + ".summary.tsv", "and", prefix + ".top5.tsv")


if __name__ == "__main__":
    main()
