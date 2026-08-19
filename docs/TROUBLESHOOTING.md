# 已知陷阱（Troubleshooting）— 本轮实战踩过的坑

## 1. PC 解释度（PC%）别算错

- `plot_panel_pca.py` / `plot_smartpca_evec.py` / `plot_smartpca_evec_groups.py` 已加**截断告警**：
  若 eigenval/eval 文件行数 ≤ 画的 PC 数，说明分母被截断（只算了前 N 个特征值），PC% 会虚高。
- **正确做法**：
  - smartpca：`.eval` 是**完整谱**（行数 ≈ 样本数−1，如 718 样本 → 718 行，`sum = n_samples - 1`），直接用；
  - plink2：`plink2 --pca` **不要带数字**（不带 `--pca 10`），让它输出全谱；带数字只有前 N 个，会截断。

## 2. metaDMG 前置要求

- **BAM 必须按 readname 排序**（`samtools sort -n`），坐标排序直接报 `Input alignment file is not sorted`。
- **read 长度必须 < 200 bp**：metaDMG 有硬断言 `l_qseq < 200`，长 read 直接 `Assertion ... Aborted`。
  滤掉：`samtools view -h -e 'qlen < 200' -o out.bam in.bam`（保持 readname 排序）。
- 单一参考（irgsp）**不需要 lca/dfit/aggregate**，`getdamage` 就够了（输出的 `pos` 表直接给出
  5' C>T、3' G>A 频率）。

## 3. FASTQ 合并：别把 .gz 二进制拼进文本流

`cat a.fq b.fastq.gz c.fastq.gz` 会把 gzip 二进制混进文本，下游 awk/pigz/bowtie2 全崩。
**先解压再合并**：

```bash
{ cat a.fq; gzip -dc b.fastq.gz; gzip -dc c.fastq.gz; } \
| LC_ALL=C awk 'NR%4==1{h=$0} NR%4==2{s=$0} NR%4==3{p=$0} NR%4==0{if(length(s)==length($0)){print h"\n"s"\n"p"\n"$0}}' \
| pigz -p 8 -c > out.fastq.gz
```

（顺带滤掉 seq/qual 长度不匹配的坏记录，避免 bowtie2 崩。）

## 4. sbatch 里 conda 激活

非交互 batch 里 `source activate <env>` 会报 `source: activate: file not found`。用完整路径：

```bash
source /home/usr/yinmt/.local/mamba/etc/profile.d/conda.sh
conda activate snakemake
export PATH=/home/usr/yinmt/software/EIG/bin:$PATH   # smartpca/convertf
```

或直接给工具绝对路径（如 `/home/apps/bowtie2/bowtie2`、`/home/apps/samtools-1.17/samtools`），不依赖 module。

## 5. reads/覆盖统计：样本名去后缀

`basename xxx.dedup.bam` 后记得剥掉 `.dedup` / `.besthit_oryza.irgsp`，否则和 .evec/.calls 的裸 ID
对不上，master 表里 reads/覆盖全变 NA。`tools/build_master_table.py` 内部已做 `norm()` 自动去后缀。

## 6. 覆盖率统计口径

- `raw_panel_covered` / `maf_covered` = 720 面板位点覆盖数（来自 `19_survey_ancient_coverage.py`）。
- 「irgsp 基因组覆盖位点数」用 `samtools depth -a bam | awk '$3>0{c++}END{print c}'`（按参考位置数，
  大 BAM 慢）。

## 7. 顺序核对（merge 前必做）

`.calls.txt` 顺序必须等于参考 `.snp` 顺序，否则静默错位：

```bash
diff <(cut -f2 <sample>.bim) <(awk '{print $1}' reference.snp) && echo ORDER_OK || echo ORDER_MISMATCH
```

## 8. 小群体（如 RAY n=9）的「最近群体」不可靠

质心样本太少会「吸」样本。用 `tools/distance_ranking.py` 的 Bootstrap（`nearest_boot_pct`）
和 `gap_d2d1` 判断归属是否稳：`boot_pct` 低或 `gap` 极小 → 归属存疑。
