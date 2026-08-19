#!/bin/bash
# Call one ancient sample at a marker set via samtools mpileup | pileupCaller --randomHaploid.
# REF/ALT are taken as BIM A2=REF, A1=ALT (the marker bfile must already be irgsp-oriented).
set -euo pipefail

usage() {
  cat <<EOF
usage: $0 --bam BAM --sample SAMPLE --bfile MARKER_PLINK --ref-fasta FASTA \
          --mapq N --baseq N --seed N --out-dir DIR --label LABEL \
          [--snp FILE --sites-bed FILE]
  --snp/--sites-bed: use a precomputed pileupCaller .snp + mpileup sites.bed
          (built once from the same bfile) instead of regenerating per call.
          Must be given together; otherwise both are derived from --bfile.
EOF
}

BAM=""; SAMPLE=""; BFILE=""; REF=""; MAPQ=""; BASEQ=""; SEED=""; OUT=""; LABEL=""
SNP_SRC=""; SITES_SRC=""
PILEUP_CALLER="${PILEUP_CALLER:-pileupCaller}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --bfile) BFILE="$2"; shift 2 ;;
    --ref-fasta) REF="$2"; shift 2 ;;
    --mapq) MAPQ="$2"; shift 2 ;;
    --baseq) BASEQ="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --out-dir) OUT="$2"; shift 2 ;;
    --label) LABEL="$2"; shift 2 ;;
    --snp) SNP_SRC="$2"; shift 2 ;;
    --sites-bed) SITES_SRC="$2"; shift 2 ;;
    *) usage; exit 2 ;;
  esac
done
for v in BAM SAMPLE BFILE REF MAPQ BASEQ SEED OUT LABEL; do
  [[ -n "${!v}" ]] || { echo "missing --${v}"; usage; exit 2; }
done

mkdir -p "$OUT"
PREFIX="$OUT/$LABEL"

# Fail fast if the input bfile is missing/empty or the working directory is wrong.
# A silently-empty .snp/.sites.bed here is exactly what produces all-9 projections.
[[ -s "${BFILE}.bim" ]] || { echo "FATAL: ${BFILE}.bim not found (pwd=$(pwd))" >&2; exit 1; }
[[ -s "${BFILE}.bed" ]] || { echo "FATAL: ${BFILE}.bed not found (pwd=$(pwd))" >&2; exit 1; }

# BIM columns: CHR SNP_ID CM BP A1 A2
# pileupCaller .snp: SNP_ID CHR CM BP REF ALT  (REF=A2, ALT=A1)
# normalize CHR "1" -> "chr01" to match the BAM/FASTA contig names.
# If --snp/--sites-bed were passed, reuse the precomputed shared files instead
# of re-deriving a 6.7M-line .snp/.sites.bed per sample (376x redundant).
if [[ -n "$SNP_SRC" || -n "$SITES_SRC" ]]; then
  [[ -n "$SNP_SRC" && -n "$SITES_SRC" ]] || { echo "FATAL: --snp and --sites-bed must be given together" >&2; exit 2; }
  echo "using precomputed snp/sites.bed: $SNP_SRC / $SITES_SRC"
  SNP_SRC="$(cd "$(dirname "$SNP_SRC")" && pwd)/$(basename "$SNP_SRC")"
  SITES_SRC="$(cd "$(dirname "$SITES_SRC")" && pwd)/$(basename "$SITES_SRC")"
else
  SNP_SRC="$PREFIX.snp"
  SITES_SRC="$PREFIX.sites.bed"
  awk 'BEGIN{OFS="\t"} {
    n=$1+0; chr=sprintf("chr%02d", n);
    print $2, chr, $3, $4, $6, $5
  }' "${BFILE}.bim" | sort -k2,2V -k4,4n > "$SNP_SRC"

  # sites.bed (0-based, half-open), sorted to match mpileup -l
  awk 'BEGIN{OFS="\t"} {print $2, $4-1, $4, $1}' "$SNP_SRC" > "$SITES_SRC"
fi

[[ -s "$SNP_SRC" ]] || { echo "FATAL: $SNP_SRC is empty or missing" >&2; exit 1; }
[[ -s "$SITES_SRC" ]] || { echo "FATAL: $SITES_SRC is empty or missing" >&2; exit 1; }

echo "=== $LABEL mapq=$MAPQ baseq=$BASEQ seed=$SEED ==="
samtools mpileup -R -B -q "$MAPQ" -Q "$BASEQ" -l "$SITES_SRC" -f "$REF" "$BAM" \
  | "$PILEUP_CALLER" --randomHaploid --seed "$SEED" --sampleNames "$SAMPLE" --samplePopName Rice \
      -f "$SNP_SRC" -p "$PREFIX" \
  2> "$PREFIX.pileupcaller.stderr"

cat "$PREFIX.pileupcaller.stderr"
echo "output: $PREFIX.bed/.bim/.fam"
