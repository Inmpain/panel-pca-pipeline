#!/usr/bin/env bash
# Submit a SLURM array job that runs one ancient sample per task through the
# validated pileupcaller_shared_call.sh + pileupcaller_plink_to_calls.py chain,
# sharing ONE precomputed marker .snp / .sites.bed across all tasks so the
# 6.7M-line marker files are derived from the bfile only once, not per sample.
#
# Dual mode (dispatch on SLURM_ARRAY_TASK_ID):
#   submitter (no SLURM_ARRAY_TASK_ID): validate inputs; build OUT/shared.snp +
#     OUT/shared.sites.bed from the A2=irgsp-locked bfile once; write the task
#     script to OUT/_slurm/; sbatch --array=1-N where N = valid samples.
#   task (SLURM_ARRAY_TASK_ID set): sample = line $SLURM_ARRAY_TASK_ID of the
#     sample list; run shared_call (with --snp/--sites-bed) + plink_to_calls
#     into OUT/q<mapq>/SAMPLE.calls.txt.
#
# Task env needs: samtools, python3, plink2, pileupCaller. The sbatch job
# inherits the submitter's environment by default (--export=ALL), so run the
# submitter from a shell where the snakemake conda env is already active.
set -euo pipefail

usage() {
  cat <<EOF
usage: $0 --sample-list FILE --bfile MARKER_PLINK --bam-dir DIR --ref-fasta FASTA \\
        --out-dir DIR [options]

  --sample-list   one ancient sample ID per line (no .bam suffix)
  --bfile         A2=irgsp-locked marker PLINK bfile (full 6.7M panel)
  --bam-dir       dir of SAMPLE<--bam-suffix>.bam
  --bam-suffix    BAM suffix (default .dedup.bam)
  --ref-fasta     reference FASTA (irgsp.fa)
  --out-dir       calls land in OUT/q<mapq>/; shared marker files at OUT/
  --mapq          minMapQ (default 25)
  --baseq         minBaseQ (default 30)
  --seed          pileupCaller seed (default 0)
  --pileup-caller path to pileupCaller binary (default \$PILEUP_CALLER or pileupCaller)
  --partition     SLURM partition (default comp)
  --exclude       SLURM --exclude nodes (default node05,node06)
  --mem           per-task memory (default 16G)
  --time          per-task wall time (default 12:00:00)
  --cpus          per-task cpus (default 1)
  --throttle      optional array throttle N -> --array=1-Nvalid%N
  --job-name      sbatch job name (default pileupcaller_q<mapq>)
  --task-prelude  literal shell lines inserted at top of each task (e.g. conda
                  activate); no shell expansion is applied to this string
  --no-submit     write the task script but do not sbatch (dry run)
EOF
}

SAMPLE_LIST=""; BFILE=""; BAMDIR=""; REF=""; OUT=""
BAM_SUFFIX=".dedup.bam"; MAPQ=25; BASEQ=30; SEED=0
PILEUP_CALLER="${PILEUP_CALLER:-pileupCaller}"
PARTITION="comp"; EXCLUDE="node05,node06"; MEM="16G"; TIME="12:00:00"; CPUS=1
THROTTLE=""; JOB_NAME=""; TASK_PRELUDE=""; DO_SUBMIT=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-list) SAMPLE_LIST="$2"; shift 2 ;;
    --bfile) BFILE="$2"; shift 2 ;;
    --bam-dir) BAMDIR="$2"; shift 2 ;;
    --bam-suffix) BAM_SUFFIX="$2"; shift 2 ;;
    --ref-fasta) REF="$2"; shift 2 ;;
    --out-dir) OUT="$2"; shift 2 ;;
    --mapq) MAPQ="$2"; shift 2 ;;
    --baseq) BASEQ="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --pileup-caller) PILEUP_CALLER="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --exclude) EXCLUDE="$2"; shift 2 ;;
    --mem) MEM="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --throttle) THROTTLE="$2"; shift 2 ;;
    --job-name) JOB_NAME="$2"; shift 2 ;;
    --task-prelude) TASK_PRELUDE="$2"; shift 2 ;;
    --no-submit) DO_SUBMIT=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done
for v in SAMPLE_LIST BFILE BAMDIR REF OUT; do
  [[ -n "${!v}" ]] || { echo "missing --${v}" >&2; usage; exit 2; }
done

# Resolve everything to absolute paths so the task script works regardless of
# the sbatch working directory. Expand a leading "~" on PILEUP_CALLER.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
abspath() {
  local p=$1
  case "$p" in
    /*) echo "$p" ;;
    *) echo "$(cd "$(dirname "$p")" && pwd)/$(basename "$p")" ;;
  esac
}
OUT="$(abspath "$OUT")"
BAMDIR="$(abspath "$BAMDIR")"
BFILE="$(abspath "$BFILE")"
REF="$(abspath "$REF")"
SAMPLE_LIST="$(abspath "$SAMPLE_LIST")"
if command -v "$PILEUP_CALLER" >/dev/null 2>&1; then
  PILEUP_CALLER="$(command -v "$PILEUP_CALLER")"
else
  PILEUP_CALLER="$(abspath "${PILEUP_CALLER/#\~/$HOME}")"
fi

[[ -s "$BFILE.bim" ]] || { echo "FATAL: ${BFILE}.bim not found" >&2; exit 1; }
[[ -s "$BFILE.bed" ]] || { echo "FATAL: ${BFILE}.bed not found" >&2; exit 1; }
[[ -s "$SAMPLE_LIST" ]] || { echo "FATAL: sample list empty/missing: $SAMPLE_LIST" >&2; exit 1; }
[[ -d "$BAMDIR" ]] || { echo "FATAL: bam dir missing: $BAMDIR" >&2; exit 1; }
[[ -s "$REF" ]] || { echo "FATAL: ref fasta missing: $REF" >&2; exit 1; }

# ---- 1. build the shared pileupCaller marker files once (same awk as
#         pileupcaller_shared_call.sh) ----
mkdir -p "$OUT"
SNP="$OUT/shared.snp"
SITES="$OUT/shared.sites.bed"
echo "building shared marker files from $BFILE.bim ..."
awk 'BEGIN{OFS="\t"} { n=$1+0; chr=sprintf("chr%02d", n); print $2, chr, $3, $4, $6, $5 }' \
  "$BFILE.bim" | sort -k2,2V -k4,4n > "$SNP"
awk 'BEGIN{OFS="\t"} { print $2, $4-1, $4, $1 }' "$SNP" > "$SITES"
[[ -s "$SNP" ]] || { echo "FATAL: shared.snp build failed" >&2; exit 1; }
[[ -s "$SITES" ]] || { echo "FATAL: shared.sites.bed build failed" >&2; exit 1; }
echo "  markers: $(wc -l < "$SNP")"
echo "  files:   $SNP / $SITES"

# ---- 2. filter the sample list to those with an existing BAM ----
SLURM_DIR="$OUT/_slurm"
mkdir -p "$SLURM_DIR"
VALID="$SLURM_DIR/samples.valid.txt"
: > "$VALID"
n=0; miss=0
while read -r S; do
  [[ -z "$S" ]] && continue
  if [[ -s "$BAMDIR/$S$BAM_SUFFIX" ]]; then
    echo "$S" >> "$VALID"; n=$((n+1))
  else
    echo "WARN: no BAM for $S -> excluded" >&2; miss=$((miss+1))
  fi
done < "$SAMPLE_LIST"
[[ $n -gt 0 ]] || { echo "FATAL: 0 valid samples in $SAMPLE_LIST" >&2; exit 1; }
echo "valid samples: $n (excluded: $miss)"

# ---- 3. write the task script: baked config via printf %q, runtime logic in a
#         fully quoted heredoc so the submitter never re-expands runtime vars ----
TASK="$SLURM_DIR/pileupcaller_array.sh"
{
  printf '#!/usr/bin/env bash\nset -euo pipefail\n'
  printf '%s\n' "$TASK_PRELUDE"
  printf 'SCRIPT_DIR=%q\n' "$SCRIPT_DIR"
  printf 'SAMPLE_LIST=%q\n' "$VALID"
  printf 'BAMDIR=%q\n' "$BAMDIR"
  printf 'BAM_SUFFIX=%q\n' "$BAM_SUFFIX"
  printf 'BFILE=%q\n' "$BFILE"
  printf 'REF=%q\n' "$REF"
  printf 'OUT=%q\n' "$OUT"
  printf 'MAPQ=%s\n' "$MAPQ"
  printf 'BASEQ=%s\n' "$BASEQ"
  printf 'SEED=%s\n' "$SEED"
  printf 'SNP=%q\n' "$SNP"
  printf 'SITES=%q\n' "$SITES"
  printf 'export PILEUP_CALLER=%q\n' "$PILEUP_CALLER"
} > "$TASK"
cat >> "$TASK" <<'TASKEOF'
module load samtools 2>/dev/null || module load samtools/ 2>/dev/null || true
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
[[ -n "$S" ]] || { echo "empty sample at index ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
BAM="$BAMDIR/$S$BAM_SUFFIX"
[[ -s "$BAM" ]] || { echo "SKIP $S (no BAM)" >&2; exit 0; }
QOUT="$OUT/q$MAPQ"
mkdir -p "$QOUT"
if [[ -s "$QOUT/$S.calls.txt" ]]; then
  echo "DONE $S (calls.txt exists)"; exit 0
fi
"$SCRIPT_DIR/pileupcaller_shared_call.sh" --bam "$BAM" --sample "$S" --bfile "$BFILE" \
  --ref-fasta "$REF" --mapq "$MAPQ" --baseq "$BASEQ" --seed "$SEED" \
  --snp "$SNP" --sites-bed "$SITES" --out-dir "$QOUT" --label "$S"
python3 "$SCRIPT_DIR/pileupcaller_plink_to_calls.py" --bfile "$QOUT/$S" --out "$QOUT/$S"
echo "OK $S"
TASKEOF
chmod +x "$TASK"

# ---- 4. submit ----
LOG_DIR="$OUT/_slurm/logs"
mkdir -p "$LOG_DIR"
ARRAY_SPEC="1-$n"
[[ -n "$THROTTLE" ]] && ARRAY_SPEC="${ARRAY_SPEC}%${THROTTLE}"
JOB_NAME="${JOB_NAME:-pileupcaller_q${MAPQ}}"
CMD=(sbatch --parsable --job-name="$JOB_NAME"
     --cpus-per-task "$CPUS" --mem "$MEM" --time "$TIME"
     --partition "$PARTITION" --exclude "$EXCLUDE"
     --output "$LOG_DIR/${JOB_NAME}_%A_%a.out"
     --error  "$LOG_DIR/${JOB_NAME}_%A_%a.err"
     --array="$ARRAY_SPEC" "$TASK")
if [[ $DO_SUBMIT -eq 1 ]]; then
  echo "submit: ${CMD[*]}"
  JOBID=$("${CMD[@]}")
  echo "submitted array $JOBID ($n tasks) -> $OUT/q$MAPQ/"
  echo "monitor: squeue -u \$USER"
else
  echo "[dry-run] task script: $TASK"
  echo "[dry-run] would run: ${CMD[*]}"
fi
