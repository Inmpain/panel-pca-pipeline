# 流程（FLOW）— v1 共享轴投影

从原始面板 + 古 BAM 到投影图，共 12 个阶段 S0–S11。所有路径/参数读 `config/config.yaml`
（本文用 `$CFG` 指代），下面的 shell 变量在每阶段开头从 config 手动填。

约定：
- `$PANEL` = 现代面板名（720/3k/...）
- `$RES` = 本流程结果根目录
- 每阶段结束后，检查输出 manifest 里的数字，再进下一阶段。

---

## S0 方向体检

确认面板 REF/ALT 与参考基因组的对齐关系（整体反标 / 混合方向 / 一致）。

```bash
python3 scripts/23_validate_snp_ref_against_fasta.py \
  --snp   <panel.snp> \
  --fasta <ref.fa> \
  --contig-format chr%02d \
  --out   $RES/panel.ref_vs_fasta.report.tsv
```

三种结论：`PASS` / `systematic_ref_alt_swap`（整体反标，可统一处理）/
`inconsistent_requires_manual_review`（混合方向，需翻链清单，见源仓库 make_flip_list.py）。

## S1 转 PLINK + 锁 A2=ref

```bash
# EIGENSTRAT -> PLINK（02 要求 .snp/.ind/.eigenstratgeno 同目录同前缀）
bash scripts/02_convert_eigenstrat_for_plink.sh \
  --dir <staging> --prefix $PANEL --out-dir $RES/plink

# 生成 SNP_ID<TAB>ref_base 参考等位基因清单
python3 scripts/make_irgsp_ref_list.py <panel.snp> <ref.fa> $RES/ref.alleles.txt

# 锁 A2 = 参考碱基（每次写新 bed 后都要重锁）
plink --bfile $RES/plink/$PANEL.plink --a2-allele $RES/ref.alleles.txt 2 1 \
  --keep-allele-order --make-bed --out $RES/plink/$PANEL.locked
```

## S2 参考样本集（建轴者 keep-list）

```bash
python3 scripts/06_build_reference_sample_set.py \
  --config $CFG --panel B --label $PANEL \
  --ind-file <panel.ind> --fam-file $RES/plink/$PANEL.locked.fam \
  --out-dir $RES/reference_sets
# -> $PANEL.reference_samples.keep
```

## S3 MAF/geno 过滤

```bash
bash scripts/07_make_fixed_markers.sh \
  --config $CFG --panel B --sensitivity primary --library-type pooled_mixed --track ALL \
  --bfile $RES/plink/$PANEL.locked --keep $RES/reference_sets/$PANEL.reference_samples.keep \
  --label $PANEL --out-dir $RES/maf_ld --stage geno_maf_only
# -> $PANEL.pooled_mixed.ALL.primary.geno_maf_filtered.{bed,bim,fam}
```

## S4 覆盖普查（古 BAM 扫面板位点）

```bash
python3 scripts/19_survey_ancient_coverage.py \
  --config $CFG --panel-snp <panel.snp> \
  --bam SAMPLE1=bam1 SAMPLE2=bam2 ... \
  --core-min-samples 1 --out-dir $RES/coverage
# -> ancient_union_sites.tsv（union 覆盖位点）+ per_sample_coverage_summary.tsv
```

## S5 覆盖候选（MAF ∩ 覆盖）

```bash
# MAF-pass 位点 ID
cut -f2 $RES/maf_ld/$PANEL.pooled_mixed.ALL.primary.geno_maf_filtered.bim | sort > /tmp/mafpass.ids
# 覆盖 union 位点 ID（19 输出 col1，跳过表头）
awk 'NR>1{print $1}' $RES/coverage/ancient_union_sites.tsv | sort > /tmp/covered.ids
comm -12 /tmp/mafpass.ids /tmp/covered.ids > $RES/covered.snplist
```

## S6 骨架（backbone 物理抽稀）

```bash
python3 scripts/08_make_5kb_thinned_markers.py \
  --config $CFG --label $PANEL \
  --geno-maf-bim $RES/maf_ld/$PANEL.pooled_mixed.ALL.primary.geno_maf_filtered.bim \
  --out-dir $RES/backbone
# -> $PANEL.paperlike_<kb>kb.fixed.snplist（window_bp 在 config.panel_B_720.paperlike_5kb）
```

## S7 合并 hybrid 面板（backbone ∪ covered）

```bash
bash scripts/union_snplists.sh \
  $RES/backbone/$PANEL.paperlike_<kb>kb.fixed.snplist \
  $RES/covered.snplist \
  $RES/hybrid.snplist
# 注意：covered 用 MAF∩coverage 还是 raw coverage，见 docs/CONFIG.md 的「hybrid 配方」小节
```

## S8 抽 hybrid bfile + 锁 A2

```bash
plink2 --bfile <panel 全位点锁 A2 bfile> --extract $RES/hybrid.snplist \
  --make-bed --out $RES/hybrid
plink --bfile $RES/hybrid --a2-allele $RES/ref.alleles.txt 2 1 \
  --keep-allele-order --make-bed --out $RES/hybrid.locked
```

> covered 若含非 MAF 位点，必须从「全位点锁 A2 bfile」抽（MAF bfile 里没有它们）。

## S9 古样本伪单倍体调用（minMapQ 扫描 + QC）

```bash
bash scripts/run_pileupcaller_mapq_matrix.sh \
  --bfile $RES/hybrid.locked \
  --samples "$SAMPLES" \
  --bam-dir <bam_dir> --ref-fasta <ref.fa> \
  --out-dir $RES/calls_matrix
# -> calls_matrix/q{0,20,25,30}/SAMPLE.calls.txt

python3 scripts/summarize_pseudohap_calls.py --calls-dir $RES/calls_matrix \
  --nmarkers "$(wc -l < $RES/hybrid.snplist)" --out $RES/ancient_qc.tsv
```

看 QC 表定 q（一般 q25 是甜点：比 q20 干净、比 q30 少丢低覆盖样本）。选中的 q 记为 `$Q`。

## S10 modern-only 诊断（口径 B，smartpca 直接出，PC% 自检）

```bash
# hybrid.locked -> EIGENSTRAT 参考三件套
bash scripts/29_convert_plink_to_eigenstrat.sh \
  --bfile $RES/hybrid.locked --out-dir $RES/eigenstrat --label $PANEL.ref

# 回填群体标签（PLINK 往返把 .ind col3 弄成 FID:IID 了）
awk 'NR==FNR{lab[$1]=$3; next} {id=$1; sub(/^[^:]*:/,"",id); print $1,$2,(id in lab?lab[id]:"NA")}' \
  <labeled.ind> $RES/eigenstrat/$PANEL.ref.ind > tmp && mv tmp $RES/eigenstrat/$PANEL.ref.ind

awk '{print $3}' $RES/eigenstrat/$PANEL.ref.ind | sort -u > $RES/eigenstrat/$PANEL.ref.poplistname.txt

bash scripts/14_run_fixed_smartpca.sh --config $CFG \
  --geno $RES/eigenstrat/$PANEL.ref.eigenstratgeno --snp $RES/eigenstrat/$PANEL.ref.snp \
  --ind  $RES/eigenstrat/$PANEL.ref.ind --poplist $RES/eigenstrat/$PANEL.ref.poplistname.txt \
  --label $PANEL.modern --out-dir $RES/eigenstrat

# PC% 自检：Σ特征值 ≈ marker 数
awk '{s+=$1} END{printf "sum=%.1f markers=%d\n", s, '"$(wc -l < $RES/hybrid.snplist)"'}' \
  $RES/eigenstrat/$PANEL.modern.eval

python3 scripts/plot_smartpca_evec.py \
  --evec $RES/eigenstrat/$PANEL.modern.evec --eval $RES/eigenstrat/$PANEL.modern.eval \
  --ind  $RES/eigenstrat/$PANEL.ref.ind --nmarkers "$(wc -l < $RES/hybrid.snplist)" \
  --title "$PANEL modern-only" --out-prefix $RES/eigenstrat/$PANEL.modern
```

肉眼确认群体结构（轴分离、无异常长尾）后再进投影。

## S11 投影（merge → lsqproject → 出图）

```bash
# 顺序核对（merge 前必做）
diff <(cut -f2 $RES/calls_matrix/q$Q/SAMPLE1.bim) \
     <(awk '{print $1}' $RES/eigenstrat/$PANEL.ref.snp) \
  && echo ORDER_OK || echo ORDER_MISMATCH

# merge：古样本 calls 拼到现代参考矩阵
python3 scripts/13_merge_ancients_fixed_panel.py \
  --reference-geno $RES/eigenstrat/$PANEL.ref.eigenstratgeno \
  --reference-ind  $RES/eigenstrat/$PANEL.ref.ind \
  --fixed-snp      $RES/eigenstrat/$PANEL.ref.snp \
  $(for S in $SAMPLES; do printf -- '--calls %s ' "$S=$RES/calls_matrix/q$Q/$S.calls.txt"; done) \
  --ancient-poplabel Ancient --label $PANEL.proj --out-dir $RES/merge

# smartpca lsqproject（输出静默，走 log，别 Ctrl-C）
bash scripts/14_run_fixed_smartpca.sh --config $CFG \
  --geno $RES/merge/$PANEL.proj.merged.eigenstratgeno --snp $RES/eigenstrat/$PANEL.ref.snp \
  --ind  $RES/merge/$PANEL.proj.merged.ind --poplist $RES/eigenstrat/$PANEL.ref.poplistname.txt \
  --label $PANEL.pca --out-dir $RES/merge

# 出图（低覆盖样本自动画三角）
python3 scripts/plot_smartpca_evec.py \
  --evec $RES/merge/$PANEL.pca.evec --eval $RES/merge/$PANEL.pca.eval \
  --ind  $RES/merge/$PANEL.proj.merged.ind --nmarkers "$(wc -l < $RES/hybrid.snplist)" \
  --title "$PANEL shared projection (pileupCaller q$Q)" \
  --out-prefix $RES/merge/$PANEL.pca
```

## 收尾自检

```bash
# 古样本是否都进了 .evec（16 个 Ancient）
awk '{print $NF}' $RES/merge/$PANEL.pca.evec | sort | uniq -c
# PC 解释度
head -10 $RES/merge/$PANEL.pca.eval
```
