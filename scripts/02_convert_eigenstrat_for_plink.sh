#!/usr/bin/env bash
# Phase 1 / Batch 1, script 02.
#
# EIGENSTRAT (.ind/.snp/.geno or .eigenstratgeno) -> PLINK PACKEDPED
# (bed/bim/fam) via EIGENSOFT convertf. Read-only against the panel's
# source directory; writes only into --out-dir. Does not filter samples
# or sites -- that is 06/07's job. This script only changes format.
#
# convertf infers input format from genotypename's file suffix, not from
# any explicit "inputformat:" parameter (confirmed in this repo's prior
# convertf work, docs/ECOTYPE_PCA_PANEL_QC_DESIGN.md section 3 / 5.6.7).
# .eigenstratgeno is the suffix that reliably means "EIGENSTRAT text" to
# convertf in this project. If the panel's geno file is named .geno instead,
# this script first sniffs whether it is plain ASCII EIGENSTRAT text; if so
# it creates a same-directory-as-out-dir symlink named *.eigenstratgeno
# pointing at the real file (never copies/mutates the source in db/), then
# feeds that symlink to convertf. If the file looks binary-packed, this
# script STOPS -- the exact packed sub-format (PACKEDANCESTRYMAP vs
# PACKEDPED) is not something convertf's suffix rule can disambiguate for
# us, and QC_DESIGN.md section 3 flagged this as unconfirmed for panel B.

set -euo pipefail

print_usage() {
  cat <<EOF
Usage: $0 --dir PANEL_DIR --prefix PREFIX --out-dir OUT_DIR [--overwrite]

  --dir       directory containing PREFIX.ind / PREFIX.snp / PREFIX.geno|.eigenstratgeno
  --prefix    file prefix, e.g. NB_final_snp / asn720.6m / civan_snp
  --out-dir   output directory for PREFIX.plink.{bed,bim,fam} and this script's log
  --overwrite allow replacing pre-existing output
EOF
}

DIR=""; PREFIX=""; OUT_DIR=""; OVERWRITE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dir) DIR="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --overwrite) OVERWRITE=1; shift ;;
    -h|--help) print_usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; print_usage; exit 1 ;;
  esac
done
if [[ -z "$DIR" || -z "$PREFIX" || -z "$OUT_DIR" ]]; then
  print_usage; exit 1
fi

mkdir -p "$OUT_DIR"
LOG="$OUT_DIR/02_convert_eigenstrat_for_plink.${PREFIX}.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== $(date -u +%FT%TZ) 02_convert_eigenstrat_for_plink.sh --prefix $PREFIX ==="

IND="$DIR/$PREFIX.ind"
SNP="$DIR/$PREFIX.snp"
GENO_ES="$DIR/$PREFIX.eigenstratgeno"
GENO_PLAIN="$DIR/$PREFIX.geno"

[[ -f "$IND" ]] || { echo "FATAL: missing $IND"; exit 2; }
[[ -f "$SNP" ]] || { echo "FATAL: missing $SNP"; exit 2; }

if [[ -f "$GENO_ES" ]]; then
  GENO_INPUT="$GENO_ES"
  echo "using existing .eigenstratgeno: $GENO_INPUT"
elif [[ -f "$GENO_PLAIN" ]]; then
  echo "no .eigenstratgeno found, sniffing $GENO_PLAIN"
  FIRSTBYTES=$(head -c 200 "$GENO_PLAIN" | tr -d '\n')
  if [[ "$FIRSTBYTES" =~ ^[0129]+$ ]]; then
    echo "sniff result: plain ASCII EIGENSTRAT genotype text (only 0/1/2/9 chars in first 200 bytes)"
    LINK="$OUT_DIR/$PREFIX.eigenstratgeno"
    if [[ -e "$LINK" && $OVERWRITE -eq 0 ]]; then
      echo "FATAL: $LINK already exists, pass --overwrite to replace"; exit 3
    fi
    ln -sf "$(readlink -f "$GENO_PLAIN")" "$LINK"
    GENO_INPUT="$LINK"
    echo "created symlink $LINK -> $GENO_PLAIN"
  else
    echo "FATAL: BLOCKED. $GENO_PLAIN does not look like plain ASCII EIGENSTRAT text"
    echo "(first 200 bytes were not all 0/1/2/9). This may be a packed binary"
    echo "EIGENSOFT format (PACKEDANCESTRYMAP or similar) -- QC_DESIGN.md section 3"
    echo "flagged this as unconfirmed for the 6.7M_720 panel. Do not guess the"
    echo "format. Run: file '$GENO_PLAIN'; head -c 100 '$GENO_PLAIN' | xxd | head"
    echo "and report the output before proceeding."
    exit 3
  fi
else
  echo "FATAL: neither $GENO_ES nor $GENO_PLAIN exists"
  exit 2
fi

BED="$OUT_DIR/$PREFIX.plink.bed"
BIM="$OUT_DIR/$PREFIX.plink.bim"
FAM="$OUT_DIR/$PREFIX.plink.fam"
PAR="$OUT_DIR/par.${PREFIX}.EIGENSTRAT.PACKEDPED"
for f in "$BED" "$BIM" "$FAM" "$PAR"; do
  if [[ -e "$f" && $OVERWRITE -eq 0 ]]; then
    echo "FATAL: $f already exists, pass --overwrite to replace"; exit 3
  fi
done

cat > "$PAR" <<EOF
genotypename: $GENO_INPUT
snpname: $SNP
indivname: $IND
outputformat: PACKEDPED
genotypeoutname: $BED
snpoutname: $BIM
indivoutname: $FAM
familynames: NO
EOF
echo "wrote $PAR:"
cat "$PAR"

echo "running convertf..."
convertf -p "$PAR"

RAW_IND_N=$(wc -l < "$IND")
RAW_SNP_N=$(wc -l < "$SNP")
OUT_FAM_N=$(wc -l < "$FAM")
OUT_BIM_N=$(wc -l < "$BIM")
echo "sanity check: input .ind lines=$RAW_IND_N vs output .fam lines=$OUT_FAM_N"
echo "sanity check: input .snp lines=$RAW_SNP_N vs output .bim lines=$OUT_BIM_N"
if [[ "$RAW_IND_N" != "$OUT_FAM_N" ]]; then
  echo "FATAL: sample count changed during conversion ($RAW_IND_N -> $OUT_FAM_N), do not proceed"
  exit 4
fi
if [[ "$RAW_SNP_N" != "$OUT_BIM_N" ]]; then
  echo "FATAL: SNP count changed during conversion ($RAW_SNP_N -> $OUT_BIM_N), do not proceed"
  exit 4
fi

echo "OK: $PREFIX converted, $OUT_FAM_N samples x $OUT_BIM_N SNPs -> $BED/$BIM/$FAM"
