#!/usr/bin/env python3
"""Survey which fixed-panel SNP positions have any qualifying read coverage
across a set of ancient BAMs, without applying the frozen calling MAPQ/BaseQ
thresholds -- this is a coverage presence check, not a genotype call.

Stage 50 postmortem (2026-08-16): the prototype's 1000 markers were the first
1000 transversion SNPs in panel file order, clustered in roughly the first
0.2Mb of chr01 -- not a genome-wide sample. That is why 0/1000 was callable
even though the BAM had thousands of mapped reads spread across all 12
chromosomes. This script produces the coverage-aware marker universe that
replaces "first N in file order": every panel site that at least one ancient
BAM actually covers (ancient_union_sites.tsv), the subset covered by at least
--core-min-samples BAMs for sensitivity analysis (ancient_core_sites.tsv),
and a per-sample coverage summary. This survey deliberately does not filter
by transversion/transition -- see 20_filter_coverage_sites_to_transversions.py
for that post-filter step, which does not require re-scanning any BAM. Wiring
the result is consumed by the registered coverage-aware Stage 50 runner.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

from fixed_projection_lib import (
    format_contig, iter_panel_snp, parse_sample_paths, refuse_existing,
    sha256_file, tally_coverage, write_tsv,
)
from lib_ecotype_v2 import load_config


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--panel-snp", required=True,
                         help="raw upstream panel .snp (e.g. civan_snp.snp), not a pre-filtered fixed marker list")
    parser.add_argument("--bam", action="append", required=True, help="repeat SAMPLE=BAM_PATH, one per ancient sample")
    parser.add_argument("--contig-format", default="chr%02d")
    parser.add_argument("--survey-mapq", type=int, default=1,
                         help="deliberately lenient; distinct from the frozen ancient.mapq=30 "
                              "used by 10_call_ancient_fixed_markers.py for actual allele calling")
    parser.add_argument("--core-min-samples", type=int, required=True,
                         help="a site enters ancient_core_sites.tsv only if at least this many BAMs cover it")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        import pysam
    except ImportError:
        print("FATAL: pysam is required", file=sys.stderr)
        return 2

    cfg = load_config(args.config)
    bam_entries = parse_sample_paths(args.bam)
    if not (1 <= args.core_min_samples <= len(bam_entries)):
        print(f"FATAL: --core-min-samples must be between 1 and {len(bam_entries)}", file=sys.stderr)
        return 3

    output_dir = Path(args.out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    union_path = output_dir / "ancient_union_sites.tsv"
    core_path = output_dir / "ancient_core_sites.tsv"
    per_sample_path = output_dir / "per_sample_coverage_summary.tsv"
    manifest_path = output_dir / "survey_manifest.json"
    excluded_flags = ("unmapped", "duplicate", "secondary", "supplementary", "qcfail")

    try:
        refuse_existing([union_path, core_path, per_sample_path, manifest_path], args.overwrite)

        panel_records = list(iter_panel_snp(args.panel_snp))
        if not panel_records:
            raise ValueError(f"{args.panel_snp}: no usable 6-column .snp rows found")
        panel_positions = defaultdict(dict)
        for record in panel_records:
            contig = format_contig(args.contig_format, record["chrom"])
            panel_positions[contig][record["pos"]] = record
        panel_contigs = set(panel_positions)

        sample_covered = {}
        sample_read_n = {}
        for sample, bam_path in bam_entries:
            covered = set()
            qualifying_read_n = 0
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                references = set(bam.references)
                if not (references & panel_contigs):
                    print(f"WARNING: {sample}: no panel contig names found in BAM references; "
                          f"panel_sites_covered will be 0 for this sample", file=sys.stderr)
                for read in bam.fetch(until_eof=True):
                    if (read.is_unmapped or read.is_duplicate or read.is_secondary
                            or read.is_supplementary or read.is_qcfail):
                        continue
                    if read.mapping_quality < args.survey_mapq:
                        continue
                    qualifying_read_n += 1
                    contig = read.reference_name
                    targets = panel_positions.get(contig)
                    if not targets:
                        continue
                    for ref_pos0 in read.get_reference_positions(full_length=False):
                        pos1 = ref_pos0 + 1
                        if pos1 in targets:
                            covered.add((contig, pos1))
            sample_covered[sample] = covered
            sample_read_n[sample] = qualifying_read_n

        tally = tally_coverage(sample_covered)
        panel_site_n = len(panel_records)

        fieldnames = ["snp_id", "chrom", "pos", "n_samples_covered", "samples_covered"]
        union_rows, core_rows = [], []
        for record in panel_records:
            contig = format_contig(args.contig_format, record["chrom"])
            samples_here = tally.get((contig, record["pos"]), set())
            if not samples_here:
                continue
            row = {
                "snp_id": record["id"], "chrom": record["chrom"], "pos": record["pos"],
                "n_samples_covered": len(samples_here),
                "samples_covered": ",".join(sorted(samples_here)),
            }
            union_rows.append(row)
            if len(samples_here) >= args.core_min_samples:
                core_rows.append(row)
        write_tsv(union_path, union_rows, fieldnames)
        write_tsv(core_path, core_rows, fieldnames)

        per_sample_rows = []
        for sample, bam_path in bam_entries:
            covered_n = len(sample_covered[sample])
            per_sample_rows.append({
                "sample": sample, "bam_path": str(bam_path),
                "qualifying_read_n": sample_read_n[sample],
                "panel_sites_total": panel_site_n, "panel_sites_covered": covered_n,
                "panel_sites_covered_fraction": f"{covered_n / panel_site_n:.10f}",
            })
        write_tsv(per_sample_path, per_sample_rows,
                  ["sample", "bam_path", "qualifying_read_n", "panel_sites_total",
                   "panel_sites_covered", "panel_sites_covered_fraction"])

        manifest = {
            "schema_version": 1, "script": "19_survey_ancient_coverage.py",
            "panel_snp": str(args.panel_snp), "panel_snp_sha256": sha256_file(args.panel_snp),
            "panel_site_n": panel_site_n, "contig_format": args.contig_format,
            "survey_mapq": args.survey_mapq,
            "survey_baseq_filter": (
                "not applied -- this is a presence-of-qualifying-read check only; "
                "final allele calling still uses the frozen ancient.mapq=30/baseq=30 "
                "from config (docs/ECOTYPE_PCA_V2_SPEC.md), unchanged by this survey"
            ),
            "excluded_read_flags": list(excluded_flags),
            "core_min_samples": args.core_min_samples, "sample_n": len(bam_entries),
            "capture_bait_bed": cfg.get("inputs", {}).get("capture_bait_bed"),
            "capture_bait_intersection_applied": False,
            "bams_sha256": {sample: sha256_file(bam_path) for sample, bam_path in bam_entries},
            "union_site_n": len(union_rows), "core_site_n": len(core_rows),
        }
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    except (OSError, ValueError) as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 3

    print(f"PASS: surveyed {len(bam_entries)} BAM(s) against {panel_site_n} panel sites; "
          f"union={len(union_rows)} core(n>={args.core_min_samples})={len(core_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
