#!/usr/bin/env bash
# Script 29: reverse of 02_convert_eigenstrat_for_plink.sh -- PLINK PACKEDPED
# (bed/bim/fam) -> EIGENSTRAT (.snp/.ind/.eigenstratgeno) via EIGENSOFT
# convertf. Exists so a plink2-produced MAF/geno-filtered bfile (e.g.
# 07_make_fixed_markers.sh --stage geno_maf_only's output) can be fed to v1's
# scripts/ecotype_pca/run_sample_panel_pca.sh, which only speaks EIGENSTRAT.
#
# plink2 is fast/memory-efficient on packed binary genotypes; a pure-Python
# MAF filter over a 29M-site x ~3000-sample text EIGENSTRAT matrix (the
# first version of this idea) would take hours -- this reuses plink2's
# already-fast MAF filtering instead of reimplementing it slowly.
set -euo pipefail

print_usage() {
  cat <<EOF
Usage: $0 --bfile PLINK_PREFIX --out-dir OUT_DIR --label LABEL [--overwrite]

  --bfile     prefix of PLINK_PREFIX.{bed,bim,fam}
  --out-dir   output directory for LABEL.{snp,ind,eigenstratgeno} and this script's log
  --label     output filename prefix
  --overwrite allow replacing pre-existing output
EOF
}

BFILE=""; OUT_DIR=""; LABEL=""; OVERWRITE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bfile) BFILE="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --label) LABEL="$2"; shift 2 ;;
    --overwrite) OVERWRITE=1; shift ;;
    -h|--help) print_usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; print_usage; exit 1 ;;
  esac
done
if [[ -z "$BFILE" || -z "$OUT_DIR" || -z "$LABEL" ]]; then
  print_usage; exit 1
fi

BED="${BFILE}.bed"; BIM="${BFILE}.bim"; FAM="${BFILE}.fam"
for f in "$BED" "$BIM" "$FAM"; do
  [[ -f "$f" ]] || { echo "FATAL: missing $f" >&2; exit 2; }
done

mkdir -p "$OUT_DIR"
LOG="$OUT_DIR/29_convert_plink_to_eigenstrat.${LABEL}.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== $(date -u +%FT%TZ) 29_convert_plink_to_eigenstrat.sh --label $LABEL ==="

GENO_OUT="$OUT_DIR/${LABEL}.eigenstratgeno"
SNP_OUT="$OUT_DIR/${LABEL}.snp"
IND_OUT="$OUT_DIR/${LABEL}.ind"
PAR="$OUT_DIR/par.${LABEL}.PACKEDPED.EIGENSTRAT"
for f in "$GENO_OUT" "$SNP_OUT" "$IND_OUT" "$PAR"; do
  if [[ -e "$f" && $OVERWRITE -eq 0 ]]; then
    echo "FATAL: $f already exists, pass --overwrite to replace" >&2; exit 3
  fi
done

cat > "$PAR" <<EOF
genotypename: $BED
snpname: $BIM
indivname: $FAM
outputformat: EIGENSTRAT
genotypeoutname: $GENO_OUT
snpoutname: $SNP_OUT
indivoutname: $IND_OUT
EOF
echo "wrote $PAR:"
cat "$PAR"

echo "running convertf..."
convertf -p "$PAR"

IN_BIM_N=$(wc -l < "$BIM")
IN_FAM_N=$(wc -l < "$FAM")
OUT_SNP_N=$(wc -l < "$SNP_OUT")
OUT_IND_N=$(wc -l < "$IND_OUT")
echo "sanity check: input .bim lines=$IN_BIM_N vs output .snp lines=$OUT_SNP_N"
echo "sanity check: input .fam lines=$IN_FAM_N vs output .ind lines=$OUT_IND_N"
if [[ "$IN_BIM_N" != "$OUT_SNP_N" ]]; then
  echo "FATAL: SNP count changed during conversion ($IN_BIM_N -> $OUT_SNP_N), do not proceed" >&2
  exit 4
fi
if [[ "$IN_FAM_N" != "$OUT_IND_N" ]]; then
  echo "FATAL: sample count changed during conversion ($IN_FAM_N -> $OUT_IND_N), do not proceed" >&2
  exit 4
fi

echo "OK: $LABEL converted, $OUT_IND_N samples x $OUT_SNP_N SNPs -> $SNP_OUT/$IND_OUT/$GENO_OUT"
