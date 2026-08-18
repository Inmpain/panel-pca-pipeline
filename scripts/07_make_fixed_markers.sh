#!/usr/bin/env bash
# Phase 1 / Batch 1 (corrected), script 07.
#
# Freezes ONE coordinate system's marker set: reference samples only ->
# (capture bait intersect, if --library-type capture) -> TV/ALL split ->
# site missingness -> MAF -> LD prune -> *.fixed.snplist +
# *.marker_manifest.tsv + md5.
#
# All numeric parameters (geno/maf/ld window/r2) are resolved from config.yaml
# via lib_ecotype_v2.resolve_marker_params -- this script never accepts a raw
# --maf / --geno / --ld-window-kb override on the command line, by design,
# so a debugging session cannot accidentally drift a frozen parameter.
#
# CORRECTED (2026-08-15, GPT review of commit 10878d7), item by item:
#
# 1. --stage {fixed|geno_maf_only} replaces the old, broken `--sensitivity
#    thinning_only` special case. The old code passed SENS="thinning_only"
#    straight into resolve_marker_params(), which only recognizes
#    primary/S1/S2/S3/S4 -- it would have raised SystemExit immediately, so
#    Panel B's paperlike_5kb route (which needs the geno/MAF-filtered, pre-LD
#    intermediate) never actually worked. --sensitivity is now always one of
#    the five real values; --stage separately controls whether the script
#    stops after geno/MAF (geno_maf_only, for 08 to consume) or continues
#    through LD pruning (fixed, the normal production path).
# 7. --panel/--library-type/--track/--sensitivity/--stage are all validated
#    against their exact enum with a `case` statement immediately after
#    argument parsing -- any typo hard-fails here with a clear message. The
#    old TV/ALL python step treated anything other than the literal string
#    'TV' as ALL, and anything other than literal 'capture' as shotgun --
#    both silent-wrong-default bugs, now impossible since the enum is
#    rejected before any of that code runs (the python step below also got
#    an explicit exhaustive if/elif/else as defense in depth).
# 8. Panel A structural checks (biallelic, chr1-12, unique SNP ID, unique
#    (chrom,pos)) are now actually run (lib_ecotype_v2.validate_panel_a_bim),
#    only for --panel A, before any filtering -- hard fail with details on
#    any violation instead of silently trusting the input .bim.
# 9. geno and MAF are now two SEPARATE plink2 --make-bed calls so the
#    manifest can report after_site_missingness and after_MAF as distinct
#    numbers (spec section 6's field list), not one combined
#    "after_site_missingness_and_MAF" figure as before. A "parameters"
#    column (compact key=value string) was added alongside the individual
#    geno/maf/ld_window_kb/ld_r2 columns to literally match the field name
#    spec section 6 lists, while keeping the individual columns for
#    convenience.
# 13. -h/--help now exits 0 (it called the same usage()-then-exit-1 path as
#     a real argument error before).
# 14. Every intermediate file this script writes (step1 extract list, the
#     geno-only bed/bim/fam, the geno+MAF bed/bim/fam, the LD prune files,
#     the final pruned bed/bim/fam) is now existence-checked against
#     --overwrite before being written, not just the two final outputs.
#
# CORRECTED (2026-08-17): --library-type gains a third value, pooled_mixed,
# alongside shotgun/capture. Our ancient BAMs are pooled_mixed_capture_plus_shotgun
# (see fixed_projection_lib.POOLED_LIBRARY_TYPE, used throughout scripts
# 09-22) -- before this change, running 07 against them required passing
# --library-type shotgun (the only enum value not requiring a bait list),
# which would have silently mislabeled pooled data as shotgun in every
# output filename and manifest row. pooled_mixed is handled identically to
# shotgun in the bait-list logic below (no filtering, BASE_LIST empty) --
# only the label is now honest.
#
# Once *.fixed.snplist exists for a given (panel, library_type, track,
# sensitivity) combination, it must never be regenerated to add/drop SNPs
# because an ancient sample was added -- pass --overwrite only when you mean
# to intentionally rebuild (e.g. after fixing an implementation bug in this
# script itself, not to chase a different SNP count).

set -euo pipefail

# PLINK2 otherwise auto-detects all host CPUs (80 on the first server test),
# ignoring the small workflow allocation. Threads affect resources, not the
# frozen statistical design. Default to one and allow the SLURM runner to
# supply its allocated CPU count through this infrastructure-only variable.
PLINK_THREADS="${RICE_PCA_PLINK_THREADS:-1}"
[[ "$PLINK_THREADS" =~ ^[1-9][0-9]*$ ]] || {
  echo "FATAL: RICE_PCA_PLINK_THREADS must be a positive integer" >&2
  exit 1
}

print_usage() {
  cat <<EOF
Usage: $0 --config CFG --panel {A|B|C} --sensitivity {primary|S1|S2|S3|S4} \\
          --library-type {shotgun|capture|pooled_mixed} --track {TV|ALL} \\
          --bfile PLINK_PREFIX --keep REFERENCE_KEEP_FILE --label LABEL \\
          --out-dir OUT_DIR [--stage {fixed|geno_maf_only}] \\
          [--capture-snp-list FILE] [--overwrite]

  --stage              fixed (default): full pipeline through LD pruning,
                        produces *.fixed.snplist + *.marker_manifest.tsv.
                        geno_maf_only: stop after site-missingness+MAF, no
                        LD pruning -- produces *.geno_maf_filtered.{bed,bim,fam}
                        + *.geno_maf_manifest.tsv for 08's paperlike_5kb route.
  --capture-snp-list   required if --library-type capture: output of
                        05_intersect_panel_baits.py ({label}.capture_compatible.snp)
EOF
}

CONFIG=""; PANEL=""; SENS=""; LIBTYPE=""; TRACK=""; BFILE=""; KEEP=""; LABEL=""; OUT_DIR=""
CAPTURE_SNP_LIST=""; OVERWRITE=0; STAGE="fixed"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --panel) PANEL="$2"; shift 2 ;;
    --sensitivity) SENS="$2"; shift 2 ;;
    --library-type) LIBTYPE="$2"; shift 2 ;;
    --track) TRACK="$2"; shift 2 ;;
    --stage) STAGE="$2"; shift 2 ;;
    --bfile) BFILE="$2"; shift 2 ;;
    --keep) KEEP="$2"; shift 2 ;;
    --label) LABEL="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --capture-snp-list) CAPTURE_SNP_LIST="$2"; shift 2 ;;
    --overwrite) OVERWRITE=1; shift ;;
    -h|--help) print_usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; print_usage; exit 1 ;;
  esac
done
for v in CONFIG PANEL SENS LIBTYPE TRACK BFILE KEEP LABEL OUT_DIR; do
  [[ -z "${!v}" ]] && { echo "missing required --${v,,}" >&2; print_usage; exit 1; }
done

# --- item 7: strict enum validation, hard fail on any typo ---
case "$PANEL" in A|B|C) ;; *) echo "FATAL: --panel must be A, B, or C, got '$PANEL'" >&2; exit 1 ;; esac
case "$LIBTYPE" in shotgun|capture|pooled_mixed) ;; *) echo "FATAL: --library-type must be shotgun, capture, or pooled_mixed, got '$LIBTYPE'" >&2; exit 1 ;; esac
case "$TRACK" in TV|ALL) ;; *) echo "FATAL: --track must be TV or ALL, got '$TRACK'" >&2; exit 1 ;; esac
case "$SENS" in primary|S1|S2|S3|S4) ;; *) echo "FATAL: --sensitivity must be one of primary/S1/S2/S3/S4, got '$SENS'" >&2; exit 1 ;; esac
case "$STAGE" in fixed|geno_maf_only) ;; *) echo "FATAL: --stage must be fixed or geno_maf_only, got '$STAGE'" >&2; exit 1 ;; esac
if [[ "$LIBTYPE" == "capture" && -z "$CAPTURE_SNP_LIST" ]]; then
  echo "FATAL: --library-type capture requires --capture-snp-list (05's output)" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
RUN_ID="${LABEL}.${LIBTYPE}.${TRACK}.${SENS}"
LOG="$OUT_DIR/07_make_fixed_markers.${RUN_ID}.${STAGE}.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== $(date -u +%FT%TZ) 07_make_fixed_markers.sh $RUN_ID stage=$STAGE ==="

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

FIXED_SNPLIST="$OUT_DIR/${RUN_ID}.fixed.snplist"
MANIFEST="$OUT_DIR/${RUN_ID}.marker_manifest.tsv"
GENO_MAF_OUT="$OUT_DIR/${RUN_ID}.geno_maf_filtered"
GENO_MAF_MANIFEST="$OUT_DIR/${RUN_ID}.geno_maf_manifest.tsv"
STEP1_EXTRACT="$OUT_DIR/${RUN_ID}.step1_extract.snplist"
GENO_ONLY_OUT="$OUT_DIR/${RUN_ID}.geno_filtered"

# --- item 14: existence-check every intermediate this run will produce ---
CANDIDATE_OUTPUTS=("$STEP1_EXTRACT" "${GENO_ONLY_OUT}.bed" "${GENO_ONLY_OUT}.bim" "${GENO_ONLY_OUT}.fam"
                   "${GENO_MAF_OUT}.bed" "${GENO_MAF_OUT}.bim" "${GENO_MAF_OUT}.fam" "$GENO_MAF_MANIFEST")
if [[ "$STAGE" == "fixed" ]]; then
  CANDIDATE_OUTPUTS+=("$OUT_DIR/${RUN_ID}.ld.prune.in" "$OUT_DIR/${RUN_ID}.ld.prune.out"
                      "$OUT_DIR/${RUN_ID}.pruned.bed" "$OUT_DIR/${RUN_ID}.pruned.bim" "$OUT_DIR/${RUN_ID}.pruned.fam"
                      "$FIXED_SNPLIST" "$MANIFEST")
fi
for f in "${CANDIDATE_OUTPUTS[@]}"; do
  if [[ -e "$f" && $OVERWRITE -eq 0 ]]; then
    echo "FATAL: $f already exists. Pass --overwrite only if you are intentionally"
    echo "rebuilding after fixing a bug in this script -- not to chase a different"
    echo "SNP count for any ancient sample."
    exit 3
  fi
done

PARAMS_JSON=$(python3 -c "
import sys, json
sys.path.insert(0, '$SCRIPT_DIR')
from lib_ecotype_v2 import load_config, resolve_marker_params
cfg = load_config('$CONFIG')
print(json.dumps(resolve_marker_params(cfg, '$PANEL', '$SENS')))
")
echo "resolved parameters from config: $PARAMS_JSON"
GENO=$(python3 -c "import json,sys; print(json.loads(sys.argv[1])['geno'])" "$PARAMS_JSON")
MAF=$(python3 -c "import json,sys; print(json.loads(sys.argv[1])['maf'])" "$PARAMS_JSON")
LD_WINDOW_KB=$(python3 -c "import json,sys; print(json.loads(sys.argv[1])['ld_window_kb'])" "$PARAMS_JSON")
LD_R2=$(python3 -c "import json,sys; print(json.loads(sys.argv[1])['ld_r2'])" "$PARAMS_JSON")
PARAMS_STR="geno=${GENO};maf=${MAF};ld_window_kb=${LD_WINDOW_KB};ld_r2=${LD_R2};track=${TRACK};library_type=${LIBTYPE};sensitivity=${SENS}"

RAW_SNPS=$(wc -l < "${BFILE}.bim")
echo "raw_snps (full panel .bim, before reference restriction): $RAW_SNPS"

# --- item 8: Panel A structural checks, before any filtering ---
if [[ "$PANEL" == "A" ]]; then
  echo "Panel A structural validation (biallelic, chr1-12, unique SNP ID, unique position)..."
  python3 -c "
import sys
sys.path.insert(0, '$SCRIPT_DIR')
from lib_ecotype_v2 import validate_panel_a_bim
problems = validate_panel_a_bim('${BFILE}.bim')
if problems:
    print('FATAL: Panel A structural checks failed:', file=sys.stderr)
    for p in problems:
        print('  - ' + p, file=sys.stderr)
    sys.exit(1)
print('Panel A structural checks: all clean (biallelic, chr1-12, unique ID, unique position)')
"
fi

# --- 1. reference-only, optional capture-bait restriction, TV/ALL split ---
if [[ "$LIBTYPE" == "capture" ]]; then
  BAIT_OVERLAP=$(wc -l < "$CAPTURE_SNP_LIST")
  echo "bait_overlap_snps: $BAIT_OVERLAP"
  BASE_LIST="$CAPTURE_SNP_LIST"
else
  BAIT_OVERLAP="NA"
  BASE_LIST=""
fi

python3 -c "
import sys
sys.path.insert(0, '$SCRIPT_DIR')
from lib_ecotype_v2 import is_transversion
base = set()
base_list = '$BASE_LIST'
if base_list:
    with open(base_list) as fh:
        base = set(l.strip() for l in fh if l.strip())
track = '$TRACK'
kept = 0
with open('${BFILE}.bim') as fh, open('$STEP1_EXTRACT', 'w') as out:
    for line in fh:
        parts = line.split()
        chrom, snpid, _cm, pos, a1, a2 = parts[:6]
        if base_list and snpid not in base:
            continue
        # item 7 defense-in-depth: exhaustive branch, no implicit ALL default.
        if track == 'TV':
            tv = is_transversion(a1, a2)
            keep_row = (tv is True)
        elif track == 'ALL':
            keep_row = True
        else:
            print(f'FATAL: unreachable track value {track!r}', file=sys.stderr)
            sys.exit(1)
        if not keep_row:
            continue
        out.write(snpid + chr(10))
        kept += 1
print(f'after_TV_ALL: {kept}', file=sys.stderr)
"
AFTER_TV_ALL=$(wc -l < "$STEP1_EXTRACT")
echo "after_TV_ALL: $AFTER_TV_ALL"
if [[ "$AFTER_TV_ALL" -eq 0 ]]; then
  echo "FATAL: 0 SNPs survive the ${LIBTYPE}/${TRACK} extract step. Stopping -- do not"
  echo "relax TV/ALL or the capture-compatible list to force a nonzero count."
  exit 3
fi

# --- 2. reference samples, site missingness ONLY (item 9: separate step from MAF) ---
plink2 --bfile "$BFILE" --keep "$KEEP" --extract "$STEP1_EXTRACT" \
  --geno "$GENO" --threads "$PLINK_THREADS" --make-bed --out "$GENO_ONLY_OUT"
AFTER_SITE_MISSINGNESS=$(wc -l < "${GENO_ONLY_OUT}.bim")
echo "after_site_missingness: $AFTER_SITE_MISSINGNESS"

# --- 3. MAF, on top of the site-missingness-filtered set (item 9: separate step) ---
plink2 --bfile "$GENO_ONLY_OUT" --maf "$MAF" --threads "$PLINK_THREADS" \
  --make-bed --out "$GENO_MAF_OUT"
AFTER_MAF=$(wc -l < "${GENO_MAF_OUT}.bim")
echo "after_MAF: $AFTER_MAF"

REF_N=$(wc -l < "$KEEP")
GENO_MAF_MD5=$(md5sum "${GENO_MAF_OUT}.bim" | cut -d' ' -f1)

{
  echo -e "panel\tlibrary_type\ttrack\tsensitivity\treference_samples_n\traw_snps\tbait_overlap_snps\tafter_TV_ALL\tafter_site_missingness\tafter_MAF\tafter_LD_or_thinning\tparameters\tgeno\tmaf\tld_window_kb\tld_r2\tmd5"
  echo -e "${PANEL}\t${LIBTYPE}\t${TRACK}\t${SENS}\t${REF_N}\t${RAW_SNPS}\t${BAIT_OVERLAP}\t${AFTER_TV_ALL}\t${AFTER_SITE_MISSINGNESS}\t${AFTER_MAF}\tNA\t${PARAMS_STR}\t${GENO}\t${MAF}\t${LD_WINDOW_KB}\t${LD_R2}\t${GENO_MAF_MD5}"
} > "$GENO_MAF_MANIFEST"
echo "wrote pre-LD intermediate: ${GENO_MAF_OUT}.{bed,bim,fam} and $GENO_MAF_MANIFEST"
echo "(for 08's paperlike_5kb route if panel B)"

if [[ "$STAGE" == "geno_maf_only" ]]; then
  echo "STAGE=geno_maf_only: stopping here as requested, no LD pruning run"
  exit 0
fi

# --- 4. LD prune -- PLINK2 native 2-arg kb syntax, per spec section 5 ---
# spec explicitly forbids the PLINK1.9-style 3-arg "100kb 10 0.2" form.
plink2 --bfile "$GENO_MAF_OUT" --indep-pairwise "${LD_WINDOW_KB}kb" "$LD_R2" \
  --threads "$PLINK_THREADS" --out "$OUT_DIR/${RUN_ID}.ld"
plink2 --bfile "$GENO_MAF_OUT" --extract "$OUT_DIR/${RUN_ID}.ld.prune.in" \
  --threads "$PLINK_THREADS" --make-bed --out "$OUT_DIR/${RUN_ID}.pruned"
AFTER_LD=$(wc -l < "$OUT_DIR/${RUN_ID}.pruned.bim")
echo "after_LD_or_thinning: $AFTER_LD"

cut -f2 "$OUT_DIR/${RUN_ID}.pruned.bim" > "$FIXED_SNPLIST"
MD5=$(md5sum "$FIXED_SNPLIST" | cut -d' ' -f1)

{
  echo -e "panel\tlibrary_type\ttrack\tsensitivity\treference_samples_n\traw_snps\tbait_overlap_snps\tafter_TV_ALL\tafter_site_missingness\tafter_MAF\tafter_LD_or_thinning\tparameters\tgeno\tmaf\tld_window_kb\tld_r2\tmd5"
  echo -e "${PANEL}\t${LIBTYPE}\t${TRACK}\t${SENS}\t${REF_N}\t${RAW_SNPS}\t${BAIT_OVERLAP}\t${AFTER_TV_ALL}\t${AFTER_SITE_MISSINGNESS}\t${AFTER_MAF}\t${AFTER_LD}\t${PARAMS_STR}\t${GENO}\t${MAF}\t${LD_WINDOW_KB}\t${LD_R2}\t${MD5}"
} > "$MANIFEST"

echo "OK: wrote $FIXED_SNPLIST (n=$AFTER_LD, md5=$MD5) and $MANIFEST"
