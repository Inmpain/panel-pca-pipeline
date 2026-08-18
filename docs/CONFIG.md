# 配置（CONFIG）参数说明

`config/config.yaml` 是唯一参数/路径来源，脚本不硬编码数值。换项目只改这个文件。

## 数值参数（冻结，别随手改）

| 参数 | 720 值 | 说明 |
|---|---|---|
| `panel_B_720.geno` | 0.20 | 位点缺失率过滤（720 缺失高，0.10 只留 0.76%，已放宽） |
| `panel_B_720.maf` | 0.01 | 次要等位基因频率下限 |
| `panel_B_720.paperlike_5kb.window_bp` | 15000 | backbone 物理抽稀窗口（bp） |
| `panel_B_720.paperlike_5kb.seed` | 20260814 | backbone 抽稀确定性 seed |
| `ancient.mapq / baseq` | 30 / 30 | mpileup `-q` / `-Q`（minMapQ 可扫 0/20/25/30） |
| `pca.num_pcs` | 10 | smartpca 输出 PC 数 |
| `pca.lsqproject` | true | 古样本被动投影 |
| `pca.numoutlieriter` | 0 | 不做 outlier 剔除 |

## 路径参数（换项目改这里）

| 键 | 内容 |
|---|---|
| `inputs.panel_B_720.dir/prefix/filtered_suffix` | 720 面板路径 |
| `inputs.ancient_bam_dir` | 古样本 BAM 目录 |
| `inputs.irgsp_reference_fasta` | 参考基因组 FASTA |
| `results_v2_root` | 结果根目录 |

## hybrid 配方（S7 的关键决策）

`hybrid = backbone ∪ covered`，其中 covered 有两种：

| 配方 | covered | 效果 |
|---|---|---|
| MAF ∩ coverage | 872（720） | 干净，但古样本 call 数少（低覆盖样本 2–7 个） |
| raw coverage | 5192（720） | 含非 MAF 位点（对现代是 singleton、对古是真覆盖），call 数涨 3–60× |

**实测结论（720）**：`backbone ∪ raw coverage`（25,866 marker）modern-only 结构不塌，
且古样本 call 数暴涨（LV7008416294 从 0 救到 42、LV7008416379 469→1522）。
**推荐用 raw coverage 配方**；非 MAF 位点在 smartpca 里被方差标准化放大，需用
modern-only 图确认结构后再定。

## minMapQ 选择

默认扫 `q0/q20/q25/q30`（`run_pileupcaller_mapq_matrix.sh --mapq`）。q25 通常是甜点：
- q20→q25 掉得少（保留信息量）；
- q25→q30 低覆盖样本 call 数对半掉（如 162→92）。

看 `ancient_qc.tsv` 的 `n_called` 矩阵定 q，别套「call rate > 50%」这类通用标准
（古 DNA 稀疏覆盖不适用）。
