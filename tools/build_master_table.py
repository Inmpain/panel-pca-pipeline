#!/usr/bin/env python3
"""Assemble the per-sample master table for a panel-PCA run.

Merges, for every ancient sample: project (derived from ID prefix), depth/age
(angkor metadata), reads stats, panel coverage, q25 call count, metaDMG damage,
and nearest modern population (centroid distance in PCA space) into one TSV.

Usage:
  python3 build_master_table.py EVEC META ANGKOR_META READS CALLS COVERAGE OUT

Inputs (all TSV):
  EVEC        smartpca .evec (sample coords + label) -> nearest-pop via centroids
  META        metaDMG summary: sample<TAB>5pC_to_T_pos0<TAB>3pG_to_A_pos0
  ANGKOR_META angkor_robot_library.txt (robot_sample_id col8 -> depth col4, age col5)
  READS       sample<TAB>total_reads<TAB>mapped_reads<TAB>mean_read_len
  CALLS       sample<TAB>q25_called
  COVERAGE    sample<TAB>raw_panel_covered<TAB>maf_covered   (ALL samples, both projects)

Sample IDs in READS/CALLS/COVERAGE are normalized by stripping common BAM
suffixes (.dedup, .besthit_oryza.irgsp) so they match the .evec IDs.
"""
import sys
from collections import Counter, defaultdict

EVEC, META, ANG, READS, CALLS, COVERAGE, OUT = sys.argv[1:8]


def norm(s):
    for suf in (".dedup", ".besthit_oryza.irgsp", ".bam", ".fastq.gz", ".fq"):
        if s.endswith(suf):
            s = s[: -len(suf)]
    return s


def read2(path, keycol=0):
    d = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and not f[0].startswith("sample"):
            d[norm(f[keycol])] = f
    return d


def read_angkor(path):
    d = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 8 and f[7] != "robot_sample_id@edna_robot_sample":
            d[f[7]] = (f[3], f[4])  # depth, age
    return d


# --- evec: centroids + ancient coords (PC1-PC2) ---
cent = defaultdict(lambda: [0.0, 0.0])
cnt = Counter()
anci = {}
for line in open(EVEC):
    f = line.split()
    if f[0].startswith("#") or len(f) < 3:
        continue
    lab = f[-1]
    pc1, pc2 = float(f[1]), float(f[2])
    if lab == "Ancient":
        anci[f[0]] = (pc1, pc2)
    else:
        cent[lab][0] += pc1
        cent[lab][1] += pc2
        cnt[lab] += 1
for lab in cent:
    cent[lab][0] /= cnt[lab]
    cent[lab][1] /= cnt[lab]


def nearest(pc1, pc2):
    best, bd = None, 1e18
    for lab, (x, y) in cent.items():
        d = ((pc1 - x) ** 2 + (pc2 - y) ** 2) ** 0.5
        if d < bd:
            bd, best = d, lab
    return best, round(bd, 5)


dmg = read2(META)
ang = read_angkor(ANG)
rd = read2(READS)
calls = read2(CALLS)
cov = read2(COVERAGE)

samples = sorted(set(list(anci) + list(calls) + list(dmg)))
with open(OUT, "w") as o:
    o.write("sample_id\tproject\tdepth_cm\tage_CE\ttotal_reads\tmapped_reads\tmean_read_len\t"
            "raw_panel_covered\tmaf_covered\tq25_called\tmetaDMG_5pCtoT\tmetaDMG_3pGtoA\t"
            "nearest_pop\tnearest_dist\tnote\n")
    for s in samples:
        proj = "nanzuo" if s.startswith("YWL1") else "angkor"
        d, a = ang.get(s, ("NA", "NA"))
        r = rd.get(s, ["NA"] * 3)
        c = cov.get(s, ["NA", "NA"])
        q = calls.get(s, ["NA"])[1] if s in calls else "NA"
        dm = dmg.get(s, ["NA", "NA"])
        np_, nd = nearest(*anci[s]) if s in anci else ("NA", "NA")
        note = "bad_library" if s == "YWL1-A3495" else ""
        o.write(f"{s}\t{proj}\t{d}\t{a}\t{r[0]}\t{r[1]}\t{r[2]}\t"
                f"{c[0]}\t{c[1]}\t{q}\t{dm[0]}\t{dm[1]}\t{np_}\t{nd}\t{note}\n")
print("wrote", OUT)
