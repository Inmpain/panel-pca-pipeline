#!/usr/bin/env python3
"""Validate a .snp file's declared REF/ALT columns against the reference
FASTA sequence.

Adapted from the same-lab Snakefile.pseudohaploid.from_panel's
validate_bim_a2_against_reference rule (2026-08-17 review): our own
pipeline had no equivalent check before this. fixed_projection_lib.iter_snp
only validates that REF/ALT are ACGT and distinct from each other -- it
never checks REF against the actual reference genome.

Distinguishes two very different failure shapes rather than treating any
mismatch as fatal:

- A CLEAN, TOTAL swap onto ALT (every checked site's FASTA base equals the
  declared ALT, never the declared REF, and zero contig/range problems) is
  the signature of the panel's REF/ALT columns being labeled backwards
  relative to this FASTA -- a real but harmless-for-PCA-math data quirk
  (see docs/decisions_log.md, 2026-08-17 entry, for why this does not
  corrupt any existing result: the modern .eigenstratgeno dosage and the
  ancient BAM-base-vs-label comparison are both self-consistent within the
  .snp file's own ref/alt label space and never cross-reference the true
  FASTA). This PASSes (exit 0), loudly logged as
  pattern=systematic_ref_alt_swap.
- Anything else -- a MIX of ref-matching and alt-matching sites, any
  genuine mismatch against both columns, or any contig-not-found /
  out-of-range position -- still hard-FATALs. That shape is what a truly
  corrupted panel or a wrong reference build would produce (scattered,
  inconsistent mismatches), and must not be silently waved through.

Writes a small summary report (not one row per SNP -- civan_snp.snp alone
has 2.36M rows).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from fixed_projection_lib import format_contig, iter_panel_snp, refuse_existing


def load_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, str] = {}
    chrom = None
    chunks: list[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if chrom is not None:
                    seqs[chrom] = "".join(chunks).upper()
                chrom = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if chrom is not None:
            seqs[chrom] = "".join(chunks).upper()
    if not seqs:
        raise ValueError(f"{path}: no sequences parsed (empty or malformed FASTA)")
    return seqs


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snp", required=True, help="raw upstream panel .snp, e.g. civan_snp.snp (6-column)")
    parser.add_argument("--fasta", required=True, help="reference FASTA the ancient BAMs were mapped against")
    parser.add_argument("--contig-format", default="chr%02d")
    parser.add_argument("--max-report", type=int, default=30)
    parser.add_argument("--out", required=True, help="small summary report path")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        refuse_existing([args.out], args.overwrite)
        seqs = load_fasta(args.fasta)
        checked = 0
        match_ref = 0
        match_alt = 0
        true_mismatches: list[str] = []
        no_contig: list[str] = []
        out_of_range: list[str] = []
        for record in iter_panel_snp(args.snp):
            checked += 1
            contig = format_contig(args.contig_format, record["chrom"])
            seq = seqs.get(contig)
            if seq is None:
                no_contig.append(f"{record['id']}: contig {contig!r} not found in FASTA")
                continue
            idx = record["pos"] - 1
            if idx < 0 or idx >= len(seq):
                out_of_range.append(f"{record['id']}: position {record['pos']} outside contig {contig} (len={len(seq)})")
                continue
            fasta_base = seq[idx]
            if fasta_base == record["ref"]:
                match_ref += 1
            elif fasta_base == record["alt"]:
                match_alt += 1
            else:
                true_mismatches.append(
                    f"{record['id']}: declared REF={record['ref']}/ALT={record['alt']} but "
                    f"FASTA {contig}:{record['pos']}={fasta_base}"
                )

        problem_free = not true_mismatches and not no_contig and not out_of_range
        if problem_free and match_alt == 0 and match_ref > 0:
            pattern = "matches_true_reference"
        elif problem_free and match_ref == 0 and match_alt > 0:
            pattern = "systematic_ref_alt_swap"
        else:
            pattern = "inconsistent_requires_manual_review"

        with open(args.out, "w") as out:
            out.write(f"checked_sites\t{checked}\n")
            out.write(f"match_ref\t{match_ref}\n")
            out.write(f"match_alt\t{match_alt}\n")
            out.write(f"true_mismatch_n\t{len(true_mismatches)}\n")
            out.write(f"no_such_contig_n\t{len(no_contig)}\n")
            out.write(f"out_of_range_n\t{len(out_of_range)}\n")
            out.write(f"pattern\t{pattern}\n")
            out.write(f"fasta\t{args.fasta}\n")
            out.write(f"snp\t{args.snp}\n")
            for line in true_mismatches[: args.max_report]:
                out.write(f"true_mismatch_example\t{line}\n")
            for line in no_contig[: args.max_report]:
                out.write(f"no_such_contig_example\t{line}\n")
            for line in out_of_range[: args.max_report]:
                out.write(f"out_of_range_example\t{line}\n")

        if pattern == "inconsistent_requires_manual_review":
            print(
                f"FATAL: {args.snp} vs {args.fasta} is inconsistent (not a clean ref-match or a clean "
                f"total swap): match_ref={match_ref} match_alt={match_alt} true_mismatch={len(true_mismatches)} "
                f"no_such_contig={len(no_contig)} out_of_range={len(out_of_range)} of {checked} checked. "
                f"See {args.out}.",
                file=sys.stderr,
            )
            return 3
        if pattern == "systematic_ref_alt_swap":
            print(
                f"WARNING: {args.snp}'s REF/ALT columns are systematically swapped relative to "
                f"{args.fasta} (all {match_alt}/{checked} sites match ALT, none match REF or anything else). "
                f"This is a confirmed data-labeling quirk, not corruption -- see docs/decisions_log.md. "
                f"PASSing.",
                file=sys.stderr,
            )
    except (OSError, ValueError) as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 3
    print(f"PASS: {checked} sites checked against {args.fasta}, pattern={pattern}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
