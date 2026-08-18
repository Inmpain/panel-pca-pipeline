#!/usr/bin/env python3
"""Append ancient calls to a frozen modern EIGENSTRAT matrix; missing stays 9."""

import argparse
import itertools
import json
import sys
from pathlib import Path

from fixed_projection_lib import (
    atomic_text_writer, parse_sample_paths, read_calls, read_ind, refuse_existing,
    sha256_file,
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-geno", required=True)
    parser.add_argument("--reference-ind", required=True)
    parser.add_argument("--fixed-snp", required=True)
    parser.add_argument("--calls", action="append", required=True, help="repeat SAMPLE=CALLS_PATH")
    parser.add_argument("--ancient-poplabel", default="Ancient")
    parser.add_argument("--label", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    geno_out = out_dir / f"{args.label}.merged.eigenstratgeno"
    ind_out = out_dir / f"{args.label}.merged.ind"
    manifest_out = out_dir / f"{args.label}.merge_manifest.json"
    try:
        refuse_existing([geno_out, ind_out, manifest_out], args.overwrite)
        modern = read_ind(args.reference_ind)
        modern_ids = {row["id"] for row in modern}
        entries = parse_sample_paths(args.calls)
        collisions = sorted(modern_ids & {sample for sample, _ in entries})
        if collisions:
            raise ValueError(f"ancient IDs collide with modern IDs: {collisions}")
        call_strings = {sample: read_calls(path) for sample, path in entries}
        fixed_n = sum(1 for line in open(args.fixed_snp) if line.strip())
        if fixed_n == 0:
            raise ValueError("fixed .snp is empty")
        for sample, calls in call_strings.items():
            if len(calls) != fixed_n:
                raise ValueError(f"{sample}: {len(calls)} calls != {fixed_n} markers")
            if "1" in calls:
                raise ValueError(f"{sample}: pseudo-haploid call file contains genotype 1")
        writer = atomic_text_writer(geno_out)
        rows_written = 0
        try:
            with open(args.reference_geno) as handle:
                for row_index, line in enumerate(handle):
                    if row_index >= fixed_n:
                        raise ValueError("reference genotype has more rows than fixed .snp")
                    genotype = line.strip()
                    if len(genotype) != len(modern) or set(genotype) - set("0129"):
                        raise ValueError(f"invalid reference genotype row {row_index + 1}")
                    writer.write(genotype + "".join(call_strings[sample][row_index] for sample, _ in entries) + "\n")
                    rows_written += 1
            if rows_written != fixed_n:
                raise ValueError(f"reference genotype has {rows_written} rows, fixed .snp has {fixed_n}")
        except BaseException:
            writer.abort()
            raise
        writer.commit()
        with open(ind_out, "w") as handle:
            for row in modern:
                handle.write(f"{row['id']}\t{row['sex']}\t{row['label']}\n")
            for sample, _ in entries:
                handle.write(f"{sample}\tU\t{args.ancient_poplabel}\n")
        manifest = {
            "schema_version": 1, "fixed_marker_n": fixed_n, "modern_sample_n": len(modern),
            "ancient_sample_n": len(entries), "ancient_ids": [sample for sample, _ in entries],
            "ancient_poplabel": args.ancient_poplabel, "missing_encoding": "9",
            "ancient_builds_axes": False, "fixed_snp_sha256": sha256_file(args.fixed_snp),
            "reference_geno_sha256": sha256_file(args.reference_geno),
        }
        manifest_out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    except (OSError, ValueError) as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 3
    print(f"PASS: appended {len(entries)} ancient sample(s); all {fixed_n} marker rows retained")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
