# panel-pca-pipeline

从「古 DNA BAM + 现代参考面板」到「共享坐标轴 PCA 投影图」的可复现分析流程。

- **输入**：现代参考面板（EIGENSTRAT/PLINK）+ 古样本 BAM + 参考基因组 FASTA。
- **输出**：覆盖度漏斗、位点 QC 表、modern-only 诊断图、古样本共享投影图。
- **方法**：modern-only backbone 骨架 + ancient-coverage 位点合成 hybrid 面板，
  `samtools mpileup | pileupCaller --randomHaploid` 伪单倍体调用，smartpca `lsqproject`
  投影。详见 `docs/FLOW.md`。

## 快速开始

```bash
# 1) 环境
source activate /path/to/snakemake-env   # plink / plink2 / python3 / matplotlib / yaml
module load samtools                      # samtools mpileup
export PILEUP_CALLER=~/software/pileupCaller-linux   # v1.5.3.1（1.6 会 segfault，别用）

# 2) 下载数据（百度网盘链接见 docs/INPUT.md 或下方 DATA_LINK）
#    参考基因组 FASTA + 现代面板 EIGENSTRAT 三件套

# 3) 改 config/config.yaml：把 paths 换成你的数据路径、samples 换成你的古样本

# 4) 按 docs/FLOW.md 顺序跑（或直接 bash run_pipeline.sh）
```

## 数据下载

现代参考面板 + 参考基因组存放于百度网盘，链接：

> **DATA_LINK**（待你上传后填写）

下载清单与大小见 `docs/INPUT.md`。

## 文档

| 文档 | 内容 |
|---|---|
| `docs/FLOW.md` | 阶段流程（S0–S11）+ 每一步的可执行命令 |
| `docs/INPUT.md` | 输入数据格式 + 百度网盘下载清单 |
| `docs/CONFIG.md` | config.yaml 参数说明 |
| `docs/FILE_MAP.md` | 脚本 → 阶段 → 输入/输出 对照表 |

## 目录结构

```
├── config/config.yaml       # 唯一参数/路径来源（angkor 示例，换项目改这里）
├── scripts/                 # 各阶段脚本（含 lib_ecotype_v2.py / fixed_projection_lib.py）
├── run_pipeline.sh          # bash 顺序 driver
└── docs/                    # FLOW / INPUT / CONFIG / FILE_MAP
```

## 已知陷阱（脚本已内建处理）

1. **顺序核对**：merge 前必须核对 `.calls.txt` 顺序 == 参考 `.snp` 顺序，否则静默错位。
2. **label 回填**：PLINK 往返会丢群体标签（`.ind` col3 变成 `FID:IID` 占位），
   需按样本 ID 回填（`docs/FLOW.md` S10 有命令）。
3. **PC% 自检**：modern-only 图用 smartpca `.eval` 直接算，Σ特征值应≈marker 数。

> 当前 v1 只含「共享轴投影」主链（S0–S11）。私有轴（per-sample 子集）暂未纳入。
> 脚本保留 rice（chr01/irgsp）命名，其它物种需改 `--contig-format` 与 FASTA。
