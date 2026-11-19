# ms-tournament

基于**多算法锦标赛机制**（Tournament）的质谱谱图鉴定引擎，用于从实验谱库（EDB）与预测谱库（ISDB）中，为未知代谢物的 MS/MS 谱筛选最可靠的候选化合物。

> 本引擎最初为内部平台 **MS-Analyst** 开发，现独立发布，便于同行复用与改进。

---

## 特点

- **多算法打分**：整合多种谱相似度算法（Modified Cosine、Entropy、Jaccard 等），避免单一指标造成的偏置
- **锦标赛决选**：将算法得分与匹配峰数合成协同得分（Synergy Score），每种算法提名 Top-N，再按三条轨道加权决出最终候选
- **两级质量窗口**：严格通道（Strict）优先匹配，失配时自动放宽数倍走兜底通道（Fallback），结果中标注本次匹配的来源
- **断点续传**：程序中断后可直接重跑同一命令，跳过已完成的 Query；需要全量重算时使用 `--no_resume`
- **多能量级比对**：可对 e0 / e1 / e2 三种碰撞能量全部打分并取最优
- **自带可视化**：每个成功鉴定的 Query 输出三张镜像对比图，可直接用于论文附图

---

## 工作原理

| 阶段 | 做什么 | 相关参数 |
|---|---|---|
| 0　全局预加载 | MS1 索引 CSV 与二级谱 JSON 一次性载入内存，并将字符串数值转换为数值类型，减少后续运算开销 | 5 个必填路径参数 |
| 1　MS1 快速初筛 | 按母离子质量容差从谱库中快速捞出候选化合物 | `--ppm_tol` |
| 2　多线程打分 | 每个 Query 作为独立任务进入线程池，对候选逐一计算二级谱相似度 | `--cpus`、`--compare_all_is_energies` |
| 3　锦标赛决选 | 算法得分 × 匹配峰数 → 协同得分；每条算法提名 Top-N；汇总提名后分三条轨道加权，决出最终候选 | `--min_peaks` |
| 4　精细对齐与落盘 | 对最终候选执行严格碎片对齐，失配则放宽兜底；生成汇总表与镜像对比图 | `--ms2_ppm_strict`、`--noise_threshold` |

---

## 安装

推荐用 conda 创建环境（`pyopenms` 在 conda 下更稳定）：

```bash
conda create -n ms-tournament python=3.10
conda activate ms-tournament
conda install -c bioconda pyopenms
pip install -r requirements.txt
```

若环境中已有可用的 `pyopenms`，直接执行 `pip install -r requirements.txt` 即可。

Python 3.10（开发环境实测版本）。

## 资源需求

引擎按「谱库全量载入内存」设计。参考规模：实验谱库约 103 万条、预测谱库约 47 万条，二者 MS2 JSON 原始体积合计约 2.7 GB，加载后进程峰值内存通常在 5 GB 以上。

- 建议 16 GB 内存以上
- CPU 核数越多，打分阶段提升越明显（`--cpus`）
- 谱库增大时按比例预留内存

---

## 快速开始

```bash
# 1) 准备谱库与待鉴定谱图
# 2) 最小可运行命令（参数含义见 docs/parameters.md）
#    demo_db 是精简子集，故将 --ppm_tol 放宽至 15，以保证示例谱图能命中候选
python ms_tournament.py \
  --query_mgf   examples/demo.mgf \
  --edbms1_csv  examples/demo_db/edb_ms1.csv \
  --isms1_csv   examples/demo_db/is_ms1.csv \
  --edbms2_json examples/demo_db/edb_ms2.json \
  --isms2_json  examples/demo_db/is_ms2.json \
  --output      results/demo_out \
  --ppm_tol     15 \
  --cpus        8
```

长任务建议挂后台，防止断网中断：

```bash
nohup python ms_tournament.py ... > run.log 2>&1 &
```

---

## 参数速查

| 参数 | 默认值 | 作用 |
|---|---|---|
| `--query_mgf` | 必填 | 待鉴定谱图（MGF） |
| `--edbms1_csv` / `--isms1_csv` | 必填 | EDB / ISDB 的 MS1 索引表 |
| `--edbms2_json` / `--isms2_json` | 必填 | EDB / ISDB 的 MS2 谱库 |
| `--ppm_tol` | 10.0 | MS1 母离子初筛容差（ppm） |
| `--ms2_ppm_strict` | 20.0 | MS2 碎片严格对齐容差（ppm），失配时自动放宽兜底 |
| `--noise_threshold` | 0.01 | 相对基峰强度的碎片噪声过滤阈值 |
| `--min_peaks` | 3 | 有效碎片匹配数下限，低于此值直接淘汰 |
| `--cpus` | 8 | 并行计算线程数 |
| `--output` | `ms_tournament_results` | 输出根目录 |
| `--no_resume` | 关闭 | 强制全量重算，不跳过已完成 Query |
| `--compare_all_is_energies` | 关闭 | 对全部碰撞能量级打分并取最优 |

完整说明见 [`docs/parameters.md`](docs/parameters.md)。

---

## 输出

`--output` 指定目录下：

**`tournament_summary.csv`** — 汇总表，每条候选一行。同一 Query 会输出三行，由 `selection_reason` 区分三条决选轨道：

| 分组 | 列 | 含义 |
|---|---|---|
| 标识 | `query_id` | 待鉴定谱图编号 |
| | `match_id` / `match_type` | 命中化合物 ID；来源库（`EDB` 实验库 / `ISDB` 预测库） |
| | `query_mz` / `target_mz` | 母离子 m/z 实测值与库中值 |
| | `name` / `formula` / `smiles` / `inchi` | 库中化合物信息，无则为 `N/A` |
| 逐算法打分 | `score__<算法>` / `npeak__<算法>` | 46 种指标各自的相似度与匹配峰数，见 [`docs/algorithms.md`](docs/algorithms.md) |
| | `energy` | 命中的能量级：`exp`（实验谱）或 `e0`/`e1`/`e2`（预测谱） |
| 综合 | `global_max_score` | 综合相似度得分，越接近 1 越好 |
| | `global_max_peaks` | 综合匹配峰数 |
| | `final_combined` | 三轨加权后的最终得分 |
| | `selection_reason` | 该行由哪条轨道选出：`Final_Score` / `Final_Peaks` / `Final_Combined` |
| | `matched_query_mz` / `matched_target_mz` | 实际对齐上的碎片 m/z，可用于结构解析 |
| | `extraction_track` | 本次匹配走的是严格通道（`Strict_20ppm`）还是兜底通道（`Fallback_100ppm`） |

**`Query_xxx/`** — 每个鉴定成功的谱图一个文件夹，含 `Final_Score`、`Final_Peaks`、`Final_Combined` 三张镜像对比图。

> 注意：`extraction_track` 中的数值为当前实现写死的标注（`20` / `100`）。若通过 `--ms2_ppm_strict` 调整了容差，该列仍显示默认值，请以参数为准。

---

## 数据可用性

完整谱库体量较大且涉及数据管理规定，不随代码发布。仓库中的 `examples/demo_db/` 仅提供**格式示例与极少量化合物**，用于跑通流程。

谱库格式说明见 [`docs/database_format.md`](docs/database_format.md)。使用自有谱库时，按该文档整理即可接入。

---

## 引用

若本引擎对你的研究有帮助，欢迎引用：

```
（待补：对应论文或软件说明）
```

---

## 许可

MIT License，见 [LICENSE](LICENSE)。
