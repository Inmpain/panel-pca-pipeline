#!/usr/bin/env python3
"""Summarize pseudo-haploid .calls.txt across a minMapQ sweep.

Reads OUT/q<N>/SAMPLE.calls.txt (one 0/2/9 per marker line, produced by
pileupcaller_plink_to_calls.py) and writes a TSV of callable-site counts and
call rate per (sample, q). Panel-agnostic: no panel-specific knowledge, so the
same script works for 720 and 3K with zero changes.
"""
import argparse
import glob
import os


def count_non9(path):
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.strip() != "9":
                n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calls-dir", required=True,
                    help="dir containing q0/ q20/ ... subdirs of SAMPLE.calls.txt")
    ap.add_argument("--nmarkers", type=int, required=True,
                    help="total marker count (hybrid snplist size) for call rate")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = []
    qdirs = sorted(glob.glob(os.path.join(args.calls_dir, "q*")))
    for qd in qdirs:
        q = os.path.basename(qd)
        for path in sorted(glob.glob(os.path.join(qd, "*.calls.txt"))):
            base = os.path.basename(path)
            sample = base[: -len(".calls.txt")]
            n_called = count_non9(path)
            rate = n_called / args.nmarkers if args.nmarkers else 0.0
            rows.append((sample, q, n_called, args.nmarkers, rate))

    rows.sort()
    with open(args.out, "w") as fh:
        fh.write("sample\tq\tn_called\tn_markers\tcall_rate\n")
        for sample, q, n_called, n_markers, rate in rows:
            fh.write(f"{sample}\t{q}\t{n_called}\t{n_markers}\t{rate:.6f}\n")
    print(f"wrote {args.out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
