#!/usr/bin/env python3
"""Intersect two or more plain newline-delimited SNP ID lists.

Used to combine 07_make_fixed_markers.sh's MAF/LD-pruned Civan
reference-only marker IDs with 19_survey_ancient_coverage.py's (or
20_filter_coverage_sites_to_transversions.py's) ancient-coverage marker
IDs. Operates on SNP ID strings only -- never touches REF/ALT/allele
columns, so it is safe regardless of any allele-orientation ambiguity a
PLINK bed/bim/fam round-trip (via convertf) might introduce; actual
genotype/allele handling downstream always goes back through the original,
already-validated civan_snp.snp file.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from fixed_projection_lib import refuse_existing


def read_ids(path: str) -> list[str]:
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snplist", action="append", required=True,
                         help="repeatable; path to a plain ID-per-line file")
    parser.add_argument("--out", required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        refuse_existing([args.out], args.overwrite)
        if len(args.snplist) < 2:
            raise ValueError("at least two --snplist arguments are required for an intersection")
        lists = [read_ids(path) for path in args.snplist]
        sets = [set(ids) for ids in lists]
        result = set.intersection(*sets)
        if not result:
            raise ValueError(f"intersection of {args.snplist} is empty")
        ordered = [snp_id for snp_id in lists[0] if snp_id in result]
        with open(args.out, "w") as out:
            out.write("\n".join(ordered) + "\n")
    except (OSError, ValueError) as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 3
    sizes = ", ".join(f"{Path(p).name}={len(l)}" for p, l in zip(args.snplist, lists))
    print(f"PASS: intersection of [{sizes}] -> {len(ordered)} IDs written to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
