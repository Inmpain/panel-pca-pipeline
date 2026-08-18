# 文件地图（FILE_MAP）— 脚本 → 阶段 → 输入/输出

| 阶段 | 脚本（scripts/） | 输入 | 输出 |
|---|---|---|---|
| S0 方向体检 | `23_validate_snp_ref_against_fasta.py` | panel `.snp` + FASTA | `*.ref_vs_fasta.report.tsv` |
| S1 转 PLINK+锁A2 | `02_convert_eigenstrat_for_plink.sh`、`make_irgsp_ref_list.py` | EIGENSTRAT 三件套 | 锁 A2 bfile、ref.alleles.txt |
| S2 参考样本集 | `06_build_reference_sample_set.py` | config + `.ind` + `.fam` | `*.reference_samples.keep` |
| S3 MAF/geno | `07_make_fixed_markers.sh` | 锁 A2 bfile + keep | `*.geno_maf_filtered.{bed,bim,fam}` |
| S4 覆盖普查 | `19_survey_ancient_coverage.py` | panel `.snp` + 古 BAM | `ancient_union_sites.tsv`、`per_sample_coverage_summary.tsv` |
| S5 覆盖候选 | `25_intersect_snplists.py`（或 comm） | MAF-pass IDs + covered IDs | covered.snplist |
| S6 骨架 | `08_make_5kb_thinned_markers.py` | geno_maf_filtered `.bim` | `*.paperlike_<kb>kb.fixed.snplist` |
| S7 合并 | `union_snplists.sh` | backbone + covered | hybrid.snplist |
| S8 抽 bfile | plink2/plink 命令 | hybrid.snplist + 全位点 bfile | hybrid.locked bfile |
| S9 古样本调用 | `run_pileupcaller_mapq_matrix.sh`、`pileupcaller_shared_call.sh`、`pileupcaller_plink_to_calls.py`、`summarize_pseudohap_calls.py` | hybrid bfile + 古 BAM | `calls_matrix/q*/SAMPLE.calls.txt`、ancient_qc.tsv |
| S10 modern 诊断 | `29_convert_plink_to_eigenstrat.sh`、`14_run_fixed_smartpca.sh`、`plot_smartpca_evec.py` | hybrid bfile | 参考 EIGENSTRAT + `.modern.png` |
| S11 投影 | `13_merge_ancients_fixed_panel.py`、`14_run_fixed_smartpca.sh`、`plot_smartpca_evec.py` | 参考 + calls | `.merged.*`、`.pca.evec/.eval`、投影图 |

## 依赖库

| 文件 | 被谁 import |
|---|---|
| `lib_ecotype_v2.py` | 06 / 07 / 08 / 27 |
| `fixed_projection_lib.py` | 13 / 19 / 25 / 29 等 |

## 工具链（版本）

| 工具 | 版本 | 备注 |
|---|---|---|
| pileupCaller | v1.5.3.1 | `~/software/pileupCaller-linux`；v1.6.0.0 segfault 禁用 |
| plink / plink2 | plink 1.90 / plink2 2.0 | snakemake 环境 |
| smartpca / convertf | EIGENSOFT | `~/software/EIG/bin/` |
| samtools | 任意 | `module load samtools` |

## 命名约定（保留自源仓库，v1 未改名）

脚本号沿用 `scripts/ecotype_pca_v2/` 的编号（02/06/07/08/13/14/19/23/25/29），
原因是脚本内部有 `lib_ecotype_v2` / `fixed_projection_lib` 的 import 依赖和交叉引用，
改名需同步改引用，留到 v2 做。各脚本职责见上表，不影响使用。
