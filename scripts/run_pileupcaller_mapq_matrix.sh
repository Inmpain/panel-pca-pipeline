#!/bin/bash
# Sweep ancient pseudo-haploid calling over minMapQ (default 0 20 25 30), keeping
# BaseQ frozen. One .calls.txt per sample per q. Panel-agnostic: the marker set is
# whatever A2=irgsp-locked PLINK bfile you pass (720 hybrid, or 3K later).
#
# Wraps pileupcaller_shared_call.sh (samtools mpileup -q MAPQ -Q BASEQ |
# pileupCaller --randomHaploid) + pileupcaller_plink_to_calls.py (PLINK -> .calls.txt).
# The q here is mpileup's -q (minMapQ); -Q (minBaseQ) is fixed at --baseq.
#
# Prereqs: env has plink2 + samtools + PILEUP_CALLER (see HANDOFF section 1), and
# --bfile is already the A2=irgsp-locked marker bfile at the hybrid marker set.
set -euo pipefail

usage() {
  cat <<EOF
usage: $0 --bfile MARKER_PLINK --samples LIST --bam-dir DIR --ref-fasta FASTA \\
          --out-dir DIR [--mapq '0 20 25 30'] [--baseq 30] [--seed 0]

  --samples   space-separated ancient sample IDs (no .bam suffix)
  --bam-dir   dir of *.besthit_oryza.irgsp.bam
  --mapq      space-separated minMapQ sweep (default "0 20 25 30")
  --baseq     fixed minBaseQ (default 30)
  --seed      stable pileupCaller --seed (default 0); same across q so q is the
              only varying factor per sample
output: OUT/q<N>/SAMPLE.calls.txt + pileupCaller PLINK + per-call stderr
EOF
}

BFILE=""; SAMPLES=""; BAMDIR=""; REF=""; OUT=""; MAPQ_LIST="0 20 25 30"; BASEQ=30; SEED=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bfile) BFILE="$2"; shift 2 ;;
    --samples) SAMPLES="$2"; shift 2 ;;
    --bam-dir) BAMDIR="$2"; shift 2 ;;
    --ref-fasta) REF="$2"; shift 2 ;;
    --out-dir) OUT="$2"; shift 2 ;;
    --mapq) MAPQ_LIST="$2"; shift 2 ;;
    --baseq) BASEQ="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    *) usage; exit 2 ;;
  esac
done
for v in BFILE SAMPLES BAMDIR REF OUT; do
  [[ -n "${!v}" ]] || { echo "missing --${v}"; usage; exit 2; }
done
[[ -s "${BFILE}.bim" ]] || { echo "FATAL: ${BFILE}.bim not found" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUT"

for q in $MAPQ_LIST; do
  for S in $SAMPLES; do
    BAM="$BAMDIR/$S.besthit_oryza.irgsp.bam"
    [[ -s "$BAM" ]] || { echo "SKIP $S q=$q (no BAM: $BAM)" >&2; continue; }
    echo "=== $S mapq=$q baseq=$BASEQ ==="
    "$SCRIPT_DIR/pileupcaller_shared_call.sh" \
      --bam "$BAM" --sample "$S" --bfile "$BFILE" --ref-fasta "$REF" \
      --mapq "$q" --baseq "$BASEQ" --seed "$SEED" \
      --out-dir "$OUT/q${q}" --label "$S" \
      || { echo "FAIL $S q=$q" >&2; continue; }
    python3 "$SCRIPT_DIR/pileupcaller_plink_to_calls.py" \
      --bfile "$OUT/q${q}/$S" --out "$OUT/q${q}/$S"
  done
done
echo "done. per-sample call counts: run summarize_pseudohap_calls.py --calls-dir $OUT --nmarkers <n> --out qc.tsv"
