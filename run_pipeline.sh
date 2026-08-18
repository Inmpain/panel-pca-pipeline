#!/usr/bin/env bash
# panel-pca-pipeline v1 — 共享轴投影 bash 顺序 runner（模板）
#
# 用法:
#   bash run_pipeline.sh --config config/config.yaml --stage S0..S11 [--samples "ID1 ID2 ..."]
#
# 这是「顺序 driver」：按 S0→S11 依次执行，每步读 config（唯一参数来源）。
# 两个决策点会停下来等人工确认：S7 的 covered 配方、S9 的 q 选择。
# 改项目时只改 config/config.yaml 的 paths + 传 --samples，不要改本脚本的数值。
#
# 环境（运行前自备）:
#   source activate <env>      # plink/plink2/python3/matplotlib/yaml
#   module load samtools
#   export PILEUP_CALLER=~/software/pileupCaller-linux
set -euo pipefail

die(){ echo "FATAL: $*" >&2; exit 1; }

# --- 参数解析 ---
CONFIG="config/config.yaml"; STAGE=""; SAMPLES=""; PANEL_LETTER="B"; PANEL_LABEL="720"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --stage)  STAGE="$2";  shift 2 ;;
    --samples) SAMPLES="$2"; shift 2 ;;
    --panel-letter) PANEL_LETTER="$2"; shift 2 ;;
    --panel-label)  PANEL_LABEL="$2";  shift 2 ;;
    *) die "unknown arg: $1" ;;
  esac
done
[[ -s "$CONFIG" ]] || die "config not found: $CONFIG"
[[ -n "$SAMPLES" ]] || die "--samples 未给（古样本 ID 列表）"

# --- 从 config 读路径/参数到 shell 变量（不做数值硬编码）---
eval "$(python3 - "$CONFIG" "$PANEL_LETTER" <<'PY'
import sys, yaml, shlex
cfg = yaml.safe_load(open(sys.argv[1]))
letter = sys.argv[2]
key = {"3K":"panel_A_3k","720":"panel_B_720","civan":"panel_C_civan"}.get(letter, "panel_B_720")
p = cfg["inputs"].get(key, {})
geno_ext = p.get("geno_ext", ".eigenstratgeno")
for k, v in {
  "PANEL_DIR": p.get("dir", ""),
  "PANEL_PREFIX": p.get("prefix", ""),
  "PANEL_SNP": f"{p.get('dir','')}/{p.get('prefix','')}.snp",
  "PANEL_IND": f"{p.get('dir','')}/{p.get('prefix','')}{p.get('filtered_suffix','')}.ind",
  "PANEL_GENO": f"{p.get('dir','')}/{p.get('prefix','')}{p.get('filtered_suffix','')}{geno_ext}",
  "BAMDIR": cfg["inputs"].get("ancient_bam_dir",""),
  "REF": cfg["inputs"].get("irgsp_reference_fasta",""),
  "RES": cfg.get("results_v2_root","results"),
}.items():
    print(f"{k}={shlex.quote(v)}")
PY
)"
[[ -s "$PANEL_SNP" ]] || die "panel .snp 不存在: $PANEL_SNP"
[[ -s "$REF" ]] || die "参考 FASTA 不存在: $REF"
mkdir -p "$RES"

SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$SELF/scripts"

# 各阶段（未实现的阶段打印 TODO 提示）
run_stage(){
  local s="$1"
  case "$s" in
    S0) python3 "$SCRIPTS/23_validate_snp_ref_against_fasta.py" \
          --snp "$PANEL_SNP" --fasta "$REF" --contig-format chr%02d \
          --out "$RES/panel.ref_vs_fasta.report.tsv" ;;
    S1) die "S1（02 转 PLINK + make_irgsp_ref_list 锁 A2）路径含 staging，请在 FLOW.md 照抄执行" ;;
    S2) die "S2（06 参考样本集）依赖 S1 的 .fam，请在 FLOW.md 照抄执行" ;;
    S3) die "S3（07 MAF/geno）依赖 S1 锁 A2 bfile，请在 FLOW.md 照抄执行" ;;
    S4) python3 "$SCRIPTS/19_survey_ancient_coverage.py" \
          --config "$CONFIG" --panel-snp "$PANEL_SNP" \
          --bam $(for s in $SAMPLES; do printf -- '--bam %s=%s ' "$s" "$BAMDIR/$s.besthit_oryza.irgsp.bam"; done) \
          --core-min-samples 1 --out-dir "$RES/coverage" ;;
    *) die "stage $s 未实现于本模板（见 FLOW.md）" ;;
  esac
}

if [[ -n "$STAGE" ]]; then
  run_stage "$STAGE"
else
  for s in S0 S4; do run_stage "$s"; done
  echo "机械阶段 S0/S4 完成。其余阶段含人工决策点，请按 docs/FLOW.md 顺序执行："
  echo "  S1 转 PLINK+锁A2 → S2 参考集 → S3 MAF/geno → S5 覆盖候选 → S6 backbone → S7 union"
  echo "  S8 抽 bfile → S9 调用+QC（定 q）→ S10 modern 诊断 → S11 投影"
fi
