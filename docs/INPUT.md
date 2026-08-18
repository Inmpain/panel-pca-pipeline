# 输入数据（INPUT）与下载清单

## 输入格式

| 数据 | 格式 | 说明 |
|---|---|---|
| 现代参考面板 | EIGENSTRAT 三件套（`.snp/.ind/.eigenstratgeno`）或 PLINK（`.bed/.bim/.fam`） | `.snp` 6 列：ID chr cm pos REF ALT；`.ind` 3 列：ID sex label；`.eigenstratgeno` 每行一个位点、每列一个样本（0/2/9） |
| 群体标签 | `.ind` 第 3 列，或单独 label 文件 | 建轴者标签（如 720 的 OrA–OrF、TRJ、IND…） |
| 古样本 BAM | 映射到同一参考基因组的 `.bam` + `.bai` | 每个古样本一个 BAM；本流程默认已 besthit/dedup 到参考基因组 |
| 参考基因组 | FASTA + `.fai` | mpileup `-f` 与 REF/ALT 对齐的锚点 |

## 百度网盘下载清单（上传后把链接填进 README 的 DATA_LINK）

> 路径为 angkor 服务器上的实际位置。大小用下面的命令自己确认，再决定上传哪些。

### 必需（跑通 720 主链）

```bash
du -sh \
  /home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa \
  /home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa.fai \
  /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.snp \
  /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.ind \
  /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.geno \
  /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.filtered.ind \
  /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.filtered.geno \
  /home/scratch/yinmt202607/db/asn720data/asn720.pop.fam
```

| 文件 | 用途 | 备注 |
|---|---|---|
| `irgsp.fa` + `.fai` | 参考基因组 | mpileup `-f`、REF/ALT 对齐 |
| `asn720.6m.snp/.ind/.geno` | 720 现代面板原始 EIGENSTRAT | `.geno` 是最大头（~6.7M×720 字符） |
| `asn720.6m.filtered.ind/.geno` | UNK 剔除 + 已标群体标签（718） | 标签 OrA–OrF 等 |
| `asn720data/asn720.pop.fam` | OrA–OrF 标签来源（按 ID 匹配） | 很小 |

### 3K 面板（后续跑 3K 用，可选）

```bash
du -sh \
  /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.snp \
  /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.ind \
  /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.eigenstratgeno \
  /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.filtered.ind \
  /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.filtered.eigenstratgeno
```

> `NB_final_snp.eigenstratgeno` ≈ 90 GB，是最大件，可按需分卷或暂不上传。

### 古样本 BAM（angkor 自有数据，是否公开自定）

```bash
du -sh /home/scratch/yinmt202607/gene/results/ecotype_pca/bam_irgsp/*.besthit_oryza.irgsp.bam
```

16 个 `*.besthit_oryza.irgsp.bam` + `.bai`。这是 angkor 项目的未发表数据，上传前请自行判断
是否公开；对「换项目复现」来说，这批 BAM 换成新项目的即可，不必共享。
