#!/usr/bin/env python3
"""Phase 1 / Batch 1 (corrected), script 06.

Writes the authoritative per-panel reference/axis-builder sample list.
Everything downstream (MAF, missingness, LD, fixed marker freeze, smartpca
poplistname) must read this file, not re-derive its own notion of
"reference" from the .ind labels directly.

Panel A (29M_3k): reference = IND ∪ AUS ∪ ARO ∪ TRJ ∪ TEJ exactly.
  HARD-FAILS if any of the 5 configured axis_labels has zero individuals
  (added in this correction pass -- the original only checked this for
  panel C, which was inconsistent: panel A's 5 labels are just as much a
  frozen precondition as panel C's 6).
  ADM is reported (project_labels) but never in the keep list.
Panel B (6.7M_720): reference = all samples in --ind-file except rows that
  are structurally malformed (empty id/label) or, if --smiss-file is given,
  samples with genotype-missing fraction == 1.0 (literally every genotype
  missing -- not ">= 0.999"). CORRECTED: the previous version had a
  --all-missing-threshold CLI flag defaulting to 0.999, which is an
  invented statistical threshold with no basis in the frozen spec ("全部
  genotype missing" means exactly that, all of it). Removed entirely --
  the check is now a structural equality (F_MISS == 1.0), not a tunable
  cutoff, so it cannot drift via a command-line typo or a "just this once"
  override. Never excluded for "population looks strange" reasons.
Panel C (civan): reference = the 6 configured domesticated labels exactly.
  HARD-FAILS if the total does not equal config's expected_axis_builder_n
  (595) -- does not proceed with a "close enough" count, per spec section 5.

CORRECTED (this pass): for every panel, hard-fails on any duplicate sample
ID in the source .ind (an ambiguous ID makes reference/non-reference
attribution meaningless), and, if --fam-file is given, verifies every kept
ID appears in the PLINK .fam EXACTLY ONCE (catches both "reference sample
silently absent from the genotype matrix" and "duplicate FAM entry" before
either could silently corrupt a --keep-based plink2 run downstream).

Read-only against the .ind/.fam files; writes only into --out-dir.
"""
import sys
from pathlib import Path
from collections import Counter

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_ecotype_v2 import (base_argparser, setup_logger, load_config,
                             read_eigenstrat_ind, write_manifest_tsv, check_output_not_present,
                             find_duplicate_ids)


def build_panel_A(cfg, rows, logger):
    axis_labels = cfg["panel_A_3k"]["axis_labels"]
    project_labels = set(cfg["panel_A_3k"].get("project_labels", []))
    axis_label_set = set(axis_labels)
    keep, other = [], Counter()
    for r in rows:
        if r["label"] in axis_label_set:
            keep.append(r)
        elif r["label"] in project_labels:
            other[f"project:{r['label']}"] += 1
        else:
            other[f"unclassified:{r['label']}"] += 1

    counts = Counter(r["label"] for r in keep)
    missing = [lb for lb in axis_labels if counts.get(lb, 0) == 0]
    if missing:
        logger.error(
            f"BLOCKED: panel A config axis_labels not present with n>0 in the .ind file: "
            f"{missing}. Distinct labels actually present: "
            f"{sorted(set(r['label'] for r in rows))}. Do not guess a corrected spelling -- "
            f"confirm the real label string on the server and update config only on "
            f"explicit instruction. (This check was missing in the original cut of this "
            f"script -- panel C had it, panel A did not, which was an inconsistency.)"
        )
        sys.exit(3)
    return keep, other


def build_panel_B(cfg, rows, logger, smiss_file):
    excluded = Counter()
    smiss = {}
    if smiss_file:
        with open(smiss_file) as fh:
            header = fh.readline().split()
            iid_idx = header.index("IID")
            fmiss_idx = header.index("F_MISS")
            for line in fh:
                parts = line.split()
                smiss[parts[iid_idx]] = float(parts[fmiss_idx])
    keep = []
    for r in rows:
        if not r["id"].strip() or not r["label"].strip():
            excluded["malformed_row"] += 1
            continue
        if smiss_file and r["id"] in smiss and smiss[r["id"]] >= 1.0:
            excluded["all_genotype_missing"] += 1
            continue
        keep.append(r)
    if not smiss_file:
        logger.warning("no --smiss-file given: 'all genotype missing' technical exclusions "
                        "were NOT checked for panel B. Run 03_audit_panel.py for panel B "
                        "first and pass its .smiss output here if any all-missing samples "
                        "need excluding (spec section 5: technical failures only, never "
                        "population-based exclusion).")
    return keep, excluded


def build_panel_C(cfg, rows, logger):
    axis_labels = set(cfg["panel_C_civan"]["axis_labels"])
    expected_n = cfg["panel_C_civan"]["expected_axis_builder_n"]
    keep, other = [], Counter()
    for r in rows:
        if r["label"] in axis_labels:
            keep.append(r)
        else:
            other[r["label"]] += 1
    counts = Counter(r["label"] for r in keep)
    missing = [lb for lb in axis_labels if counts.get(lb, 0) == 0]
    if missing:
        logger.error(
            f"BLOCKED: config axis_labels not present with n>0 in the .ind file: {missing}. "
            f"Distinct labels actually present: {sorted(set(r['label'] for r in rows))}. "
            f"Do not guess a corrected spelling -- confirm the real label string on the "
            f"server and update config.panel_C_civan.axis_labels only on explicit instruction."
        )
        sys.exit(3)
    if len(keep) != expected_n:
        logger.error(
            f"BLOCKED: axis-builder total = {len(keep)}, expected exactly {expected_n}. "
            f"Per-label breakdown: {dict(counts)}. STOP per spec section 17 acceptance "
            f"condition 'exactly 595 domesticated samples define axes' -- do not adjust "
            f"the label set or filters to force this number."
        )
        sys.exit(3)
    return keep, other


def load_fam_iid_to_fid(fam_file):
    """FID -> IID mapping is NOT guaranteed to be FID==IID: convertf with
    familynames: NO (as used by 02_convert_eigenstrat_for_plink.sh) assigns
    sequential-index FIDs, not the sample ID, in column 1 -- confirmed by
    actually running 02 against a synthetic panel during this correction
    pass (see chat report). A --keep file built as 'sampleID sampleID' would
    then silently fail to match plink2's FID+IID pairing. Returns
    {iid: [fid, fid, ...]} (list, so duplicate-IID FAM rows are detectable)."""
    iid_to_fids = {}
    with open(fam_file) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            fid, iid = parts[0], parts[1]
            iid_to_fids.setdefault(iid, []).append(fid)
    return iid_to_fids


def check_fam_correspondence(keep_ids, iid_to_fids, fam_file, logger):
    missing_from_fam = [i for i in keep_ids if len(iid_to_fids.get(i, [])) == 0]
    dup_in_fam = [i for i in keep_ids if len(iid_to_fids.get(i, [])) > 1]
    if missing_from_fam:
        sample = missing_from_fam[:10]
        logger.error(f"BLOCKED: {len(missing_from_fam)} reference sample IDs are not present "
                     f"in {fam_file} at all, e.g. {sample}. A reference sample absent from "
                     f"the genotype matrix cannot be used for MAF/missingness/LD -- this "
                     f"needs to be resolved (wrong .ind/.fam pairing? sample dropped during "
                     f"conversion?), not silently excluded.")
        sys.exit(3)
    if dup_in_fam:
        sample = dup_in_fam[:10]
        logger.error(f"BLOCKED: {len(dup_in_fam)} reference sample IDs appear MORE THAN ONCE "
                     f"in {fam_file}, e.g. {sample}. plink2 --keep cannot unambiguously select "
                     f"one of several same-ID FAM rows -- this needs to be resolved, not "
                     f"silently picking the first match.")
        sys.exit(3)
    logger.info(f"fam correspondence OK: all {len(keep_ids)} reference IDs appear exactly "
                f"once in {fam_file}")


def main():
    ap = base_argparser(__doc__)
    ap.add_argument("--panel", required=True, choices=["A", "B", "C"])
    ap.add_argument("--label", required=True)
    ap.add_argument("--ind-file", required=True)
    ap.add_argument("--fam-file", default=None,
                     help="panel's PLINK .fam (from script 02) -- if given, verifies every "
                          "kept ID corresponds to exactly one FAM row")
    ap.add_argument("--smiss-file", default=None,
                     help="panel B only: 03_audit_panel.py's .smiss output, for detecting "
                          "all-genotype-missing (F_MISS==1.0) technical exclusions")
    args = ap.parse_args()
    logger, _ = setup_logger(f"06_build_reference_sample_set.{args.label}", args.out_dir)
    cfg = load_config(args.config)

    out_dir = Path(args.out_dir)
    keep_path = out_dir / f"{args.label}.reference_samples.keep"
    manifest_path = out_dir / f"{args.label}.reference_samples_manifest.tsv"
    check_output_not_present([keep_path, manifest_path], args.overwrite, logger)

    rows = read_eigenstrat_ind(args.ind_file)
    logger.info(f"panel {args.panel} ({args.label}): {len(rows)} total rows in {args.ind_file}")

    dup_ids = find_duplicate_ids(r["id"] for r in rows)
    if dup_ids:
        sample = list(dup_ids.items())[:10]
        logger.error(f"BLOCKED: {len(dup_ids)} duplicate sample IDs in {args.ind_file}, "
                     f"e.g. {sample}. Reference/non-reference attribution is ambiguous with "
                     f"duplicate IDs -- resolve the source .ind before building a reference "
                     f"set from it, do not silently keep one occurrence.")
        sys.exit(3)

    if args.panel == "A":
        keep, other = build_panel_A(cfg, rows, logger)
    elif args.panel == "B":
        keep, other = build_panel_B(cfg, rows, logger, args.smiss_file)
    else:
        keep, other = build_panel_C(cfg, rows, logger)

    keep_ids = [r["id"] for r in keep]
    if args.fam_file:
        iid_to_fids = load_fam_iid_to_fid(args.fam_file)
        check_fam_correspondence(keep_ids, iid_to_fids, args.fam_file, logger)
        with open(keep_path, "w") as fh:
            for i in keep_ids:
                fh.write(f"{iid_to_fids[i][0]}\t{i}\n")  # real FID from FAM, not assumed FID=IID
        logger.info("keep-list FID column taken from the real .fam (not assumed equal to IID)")
    else:
        logger.warning("no --fam-file given: reference-IDs-vs-genotype-matrix correspondence "
                        "was NOT checked, and the keep-list FID column below is a GUESS "
                        "(FID=IID). This is confirmed wrong for panels converted by "
                        "02_convert_eigenstrat_for_plink.sh, whose convertf par file "
                        "(familynames: NO) assigns sequential-index FIDs, not the sample ID -- "
                        "do not use this keep-list with plink2 --keep without --fam-file.")
        with open(keep_path, "w") as fh:
            for i in keep_ids:
                fh.write(f"{i}\t{i}\n")

    counts = Counter(r["label"] for r in keep)
    manifest_rows = [{"label": lb, "n": n, "role": "reference"} for lb, n in sorted(counts.items())]
    manifest_rows += [{"label": k, "n": n, "role": "excluded_or_other"} for k, n in sorted(other.items())]
    manifest_rows.append({"label": "TOTAL_REFERENCE", "n": len(keep), "role": "reference"})
    write_manifest_tsv(manifest_path, manifest_rows, ["label", "n", "role"])

    logger.info(f"reference/axis-builder N = {len(keep)}, per-label: {dict(counts)}")
    logger.info(f"other (not reference): {dict(other)}")
    logger.info(f"wrote {keep_path} and {manifest_path}")


if __name__ == "__main__":
    main()
