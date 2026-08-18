#!/usr/bin/env python3
"""Phase 1 / Batch 1 (corrected), script 08.

Panel B "paperlike_5kb" route (spec section 5, B-paperlike): from a
site-missingness/MAF-filtered (but NOT LD-pruned) SNP set, keep at most one
SNP per non-overlapping 5000bp window, deterministically, seed=20260814.

Input must be the .geno_maf_filtered.bim produced by
07_make_fixed_markers.sh --stage geno_maf_only (i.e. geno<=0.10, maf>=0.01
already applied, no LD pruning). This script does not re-apply those
filters itself -- it only thins.

Output is named *.paperlike_5kb.fixed.snplist -- never Wang2017_exact, per
spec's explicit naming prohibition (this is a reprocessing of a different,
higher-density matrix, not a reproduction of the published analysis).

CORRECTED (2026-08-15, GPT review of commit 10878d7):

5. Window assignment was `pos // window_bp`, wrong for PLINK's 1-based
   coordinates: position 5000 (the last base of window 1..5000) would have
   landed in window 1 (5000..9999) instead of window 0. Fixed to
   `(pos-1) // window_bp` via lib_ecotype_v2.genomic_window_index, which is
   also unit-tested directly (see test_lib_ecotype_v2.py). Also reimplemented
   as a single streaming pass over the (sorted) .bim -- the previous version
   loaded every (pos, snpid) pair for the whole panel into a dict of lists
   before processing anything, which is exactly the unbounded-memory pattern
   this correction batch was told to eliminate. The streaming version holds
   only the current window's candidate SNPs in memory at once, and requires
   the input to already be sorted by chrom then position (PLINK2 --make-bed
   output always is) -- it hard-fails immediately if it ever sees position
   go backwards within a chromosome, rather than silently re-sorting (which
   would defeat the point of streaming) or silently producing a wrong result.
9. Manifest now optionally inherits panel/library_type/track/
   reference_samples_n/raw_snps/bait_overlap_snps/after_TV_ALL/
   after_site_missingness/after_MAF from 07's --stage geno_maf_only manifest
   (via --upstream-manifest), so the final row satisfies the same schema as
   07's *.marker_manifest.tsv (spec section 6) instead of an ad hoc field
   set. Without --upstream-manifest those fields are written as NA with a
   logged warning -- never fabricated.
"""
import csv
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_ecotype_v2 import (base_argparser, setup_logger, load_config, md5_file,
                             check_output_not_present, genomic_window_index)

MANIFEST_FIELDS = ["panel", "library_type", "track", "sensitivity", "reference_samples_n",
                   "raw_snps", "bait_overlap_snps", "after_TV_ALL", "after_site_missingness",
                   "after_MAF", "after_LD_or_thinning", "parameters", "md5"]


def read_upstream_manifest(path):
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    if len(rows) != 1:
        raise ValueError(f"{path}: expected exactly 1 data row, found {len(rows)}")
    return rows[0]


def stream_thin(bim_path, window_bp, seed, logger):
    """Single pass over a chrom-then-position-sorted .bim. O(candidates in one
    window) memory, not O(total SNPs). Hard-fails (not silently re-sorts) if
    the file is not actually sorted that way."""
    rng = random.Random(seed)
    kept = []
    current_chrom = None
    current_window = None
    candidates = []
    prev_pos = None
    n_in = 0

    def flush():
        if candidates:
            ordered = sorted(candidates)  # stable order before the RNG draw
            chosen = ordered[rng.randrange(len(ordered))]
            kept.append((current_chrom, current_window, chosen[1]))

    with open(bim_path) as fh:
        for lineno, line in enumerate(fh, 1):
            parts = line.split()
            chrom, snpid, pos = parts[0], parts[1], int(parts[3])
            n_in += 1
            win = genomic_window_index(pos, window_bp)

            if chrom != current_chrom:
                flush()
                current_chrom, current_window, candidates = chrom, win, []
                prev_pos = None
            elif pos < (prev_pos or 0):
                raise ValueError(
                    f"{bim_path}:{lineno}: position {pos} < previous position {prev_pos} "
                    f"within chrom {chrom} -- input is not sorted by position. Streaming "
                    f"low-memory thinning requires sorted input (PLINK2 --make-bed output "
                    f"always is); re-derive this .bim rather than have this script silently "
                    f"re-sort or produce results from unsorted data."
                )
            elif win != current_window:
                flush()
                current_window, candidates = win, []

            candidates.append((pos, snpid))
            prev_pos = pos
        flush()

    logger.info(f"streamed {n_in} input SNPs, {len(kept)} windows populated")
    return kept, n_in


def main():
    ap = base_argparser(__doc__)
    ap.add_argument("--label", required=True)
    ap.add_argument("--geno-maf-bim", required=True,
                     help="output of 07_make_fixed_markers.sh --stage geno_maf_only")
    ap.add_argument("--upstream-manifest", default=None,
                     help="07's *.geno_maf_manifest.tsv -- if given, inherits panel/"
                          "library_type/track/reference_samples_n/raw_snps/bait_overlap_snps/"
                          "after_TV_ALL/after_site_missingness/after_MAF into this script's "
                          "manifest row so it matches spec section 6's schema")
    args = ap.parse_args()
    logger, _ = setup_logger(f"08_make_5kb_thinned_markers.{args.label}", args.out_dir)
    cfg = load_config(args.config)

    window_bp = cfg["panel_B_720"]["paperlike_5kb"]["window_bp"]
    seed = cfg["panel_B_720"]["paperlike_5kb"]["seed"]
    window_kb = window_bp // 1000
    route_label = f"paperlike_{window_kb}kb"
    logger.info(f"window_bp={window_bp} seed={seed} route_label={route_label} "
                f"(from config.panel_B_720.paperlike_5kb)")

    out_dir = Path(args.out_dir)
    out_snplist = out_dir / f"{args.label}.{route_label}.fixed.snplist"
    out_manifest = out_dir / f"{args.label}.{route_label}.marker_manifest.tsv"
    check_output_not_present([out_snplist, out_manifest], args.overwrite, logger)

    try:
        kept, n_in = stream_thin(args.geno_maf_bim, window_bp, seed, logger)
    except ValueError as e:
        logger.error(f"FATAL: {e}")
        sys.exit(2)

    with open(out_snplist, "w") as fh:
        for _chrom, _win, snpid in kept:
            fh.write(snpid + "\n")
    md5 = md5_file(out_snplist)

    upstream = {}
    if args.upstream_manifest:
        upstream = read_upstream_manifest(args.upstream_manifest)
        logger.info(f"inherited upstream manifest fields from {args.upstream_manifest}: "
                    f"panel={upstream.get('panel')} library_type={upstream.get('library_type')} "
                    f"track={upstream.get('track')}")
    else:
        logger.warning("no --upstream-manifest given -- panel/library_type/track/"
                        "reference_samples_n/raw_snps/bait_overlap_snps/after_TV_ALL/"
                        "after_site_missingness/after_MAF will be written as NA, not "
                        "fabricated from this script's own (much narrower) knowledge")

    row = {
        "panel": upstream.get("panel", "NA"),
        "library_type": upstream.get("library_type", "NA"),
        "track": upstream.get("track", "NA"),
        "sensitivity": route_label,
        "reference_samples_n": upstream.get("reference_samples_n", "NA"),
        "raw_snps": upstream.get("raw_snps", "NA"),
        "bait_overlap_snps": upstream.get("bait_overlap_snps", "NA"),
        "after_TV_ALL": upstream.get("after_TV_ALL", "NA"),
        "after_site_missingness": upstream.get("after_site_missingness", "NA"),
        "after_MAF": upstream.get("after_MAF", str(n_in)),
        "after_LD_or_thinning": str(len(kept)),
        "parameters": f"window_bp={window_bp};seed={seed};input_snps={n_in}",
        "md5": md5,
    }
    with open(out_manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        w.writeheader()
        w.writerow(row)

    logger.info(f"selected {len(kept)} SNPs (1 per populated {window_bp}bp window) "
                f"from {n_in} input SNPs")
    logger.info(f"wrote {out_snplist} (md5={md5}) and {out_manifest}")


if __name__ == "__main__":
    main()
