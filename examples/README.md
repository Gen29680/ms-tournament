# examples

用于跑通流程的最小示例数据。

## demo.mgf

5 个待鉴定谱图，`Query_001` 到 `Query_005`，母离子分别为 319.2339、449.2239、462.1716、561.3603、575.3763，峰数 21 / 653 / 82 / 159 / 79。

已脱敏：去除了 MZmine 导出的 `FEATURE_ID`、`SCANS`、`RTINSECONDS`、`MERGED_STATS` 等字段，仅保留 `TITLE`、`PEPMASS`、`CHARGE`、`MSLEVEL` 与碎片峰。强度保留原始科学计数写法（如 `1.1E2`）。

## demo_db

精简谱库，每库 30 个化合物（`is_ms2.json` 实际 28 条），按上述 5 个母离子就近选取，用于让示例谱图能够命中候选。

| 文件 | 内容 |
|---|---|
| `edb_ms1.csv` | 实验谱库 MS1 索引，列：`id, pepmass, smiles` |
| `edb_ms2.json` | 实验谱库 MS2 数据 |
| `is_ms1.csv` | 预测谱库 MS1 索引，列：`id, formula, exactmass, smiles, inchi, inchikey, m+h, m+nh4, m+na` |
| `is_ms2.json` | 预测谱库 MS2 数据，含三个碰撞能量 |

> 完整谱库（实验库约 103 万条、预测库约 47 万条）不随代码发布。本目录仅用于验证流程可跑通，检索结果不代表真实鉴定结论。

因为谱库覆盖稀疏，quick start 中把 `--ppm_tol` 放宽至 15，使 5 个示例谱图都能命中。使用完整谱库时应恢复默认值 10。
