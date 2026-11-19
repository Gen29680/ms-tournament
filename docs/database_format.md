# 谱库文件格式

MS/MS Tournament 依赖两类谱库：**实验谱库（EDB）** 与 **预测谱库（ISDB）**。每类各由一对文件描述——一个 MS1 索引表（CSV，用于母离子初筛）和一个 MS2 谱库（JSON，用于碎片谱打分）。

两个文件之间通过 `id` 关联。所有文件在启动时一次性载入内存。

```
EDB:  edb_ms1.csv  (母离子索引)   +  edb_ms2.json  (碎片谱)
ISDB: is_ms1.csv   (母离子索引)   +  is_ms2.json   (碎片谱 × 3 个能量)
```

---

## 一、EDB MS1 索引表（CSV）

三列，表头固定：

| 列名 | 类型 | 说明 |
| --- | --- | --- |
| `id` | str | 化合物唯一标识，与 `edb_ms2.json` 的 key 对应，例如 `cur228946` |
| `pepmass` | float | 母离子质荷比（m/z） |
| `smiles` | str | 结构 SMILES |

示例：

```
id,pepmass,smiles
cur228946,319.238,CCCC[C@@H](C=O)NC(=O)[C@H](CC(C)C)NC(=O)OCC1=CC=CC=C1
```

## 二、EDB MS2 谱库（JSON）

顶层是 **dict**，key 为 `id`，value 为另一层 dict：

```json
{
  "cur79799": {
    "pepmass": 561.3603,
    "charge": 1,
    "ion_mode": "",
    "smiles": "CCCCC/C=C\\C/C=C\\CCCCCCCC(=O)O...",
    "compound_name": "…",
    "ms2": [[84.0806, 5.49], [201.1236, 13.39], [243.1341, 99.5], ...]
  }
}
```

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `pepmass` | float | 母离子 m/z |
| `charge` | int | 电荷数 |
| `ion_mode` | str | 离子模式，可为空字符串 |
| `smiles` | str | 结构 SMILES |
| `compound_name` | str | 化合物名 |
| `ms2` | list 或 str | 碎片峰。**两种写法代码都支持**：原生二维数组 `[[mz, 强度], ...]`，或该数组的字符串形式 `"[[84.08, 5.49], ...]"` |

碎片强度为 **0–1000 归一化**（谱图中最高峰约为 999）。

> 读取逻辑见 `my_packages/functions_new.py` 的 `get_spectra_from_einfo()`：若 `ms2` 是字符串，先用 `ast.literal_eval` 还原；是原生数组则直接使用。随后按 `max_mz = pepmass + 0.01` 清洗谱图。

## 三、ISDB MS1 索引表（CSV）

九列，表头固定：

| 列名 | 类型 | 说明 |
| --- | --- | --- |
| `id` | str | 化合物标识，例如 `CNP0381886` |
| `formula` | str | 分子式 |
| `exactmass` | float | 精确质量（中性分子） |
| `smiles` | str | 结构 SMILES |
| `inchi` | str | InChI，可为空 |
| `inchikey` | str | InChIKey |
| `m+h` | float | `[M+H]⁺` 质量 |
| `m+nh4` | float | `[M+NH₄]⁺` 质量 |
| `m+na` | float | `[M+Na]⁺` 质量 |

母离子初筛使用 `exactmass` 列。

示例：

```
id,formula,exactmass,smiles,inchi,inchikey,m+h,m+nh4,m+na
CNP0381886,C17H29N5O,319.23721,N=1C=C(N(C1)C)CN2CC3(N4CCOCC4CN(C3)C(C)C)C2,,PAKWQBNLVWCDMR-UHFFFAOYSA-N,...
```

## 四、ISDB MS2 谱库（JSON）

顶层同样是 dict，key 为 `id`，value 含三个碰撞能量下的预测碎片谱：

```json
{
  "NPA019095": {
    "smiles": "...",
    "pepmass": 575.3763,
    "energy0_ms2": "[[478.3527, 6.3], [540.36835, 8.25], [558.37891, 38.1], [576.38948, 100.0]]",
    "energy1_ms2": "[[41.03858, 10.81], [55.05423, 38.9], ...]",
    "energy2_ms2": "[[41.03858, 38.45], [55.05423, 49.9], ...]"
  }
}
```

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `smiles` | str | 结构 SMILES |
| `pepmass` | float | 母离子 m/z |
| `energy0_ms2` | **str** | 10 eV 预测碎片谱，`[[mz, 强度], ...]` 的字符串形式 |
| `energy1_ms2` | **str** | 20 eV |
| `energy2_ms2` | **str** | 40 eV |

三点必须注意：

1. **三个 energy 字段必须是字符串**，不能是原生数组。代码直接对其调用 `ast.literal_eval`，写成数组会报错。三个碰撞能量分别对应 CFM-ID 在 10 / 20 / 40 eV 的预测结果。
2. 碎片强度为 **0–100 归一化**，与 EDB 的 0–1000 不同。
3. 代码在构造预测谱图时使用 `precursor_mz = pepmass + 0.1`，这是实现中的固定偏置，自建库时无需额外处理。

> 读取逻辑见 `my_packages/functions_new.py` 的 `get_spectra_from_isinfo()`。

## 五、自建谱库的最小要求

若要用自己的数据替换 demo 库，需准备四个文件，且满足：

- **MS1 索引表**为带表头的 CSV，列名与上文一致；`id` 列的取值必须与对应 MS2 JSON 的 key 完全一致。
- **MS2 JSON** 顶层为 `{id: {...}}` 的 dict；EDB 走 `ms2` 字段，ISDB 走 `energy0_ms2` / `energy1_ms2` / `energy2_ms2` 三个字段。
- 两个库可以只提供其中一个：把另一个的参数指向同一对文件，或用极少条目构造空库即可（程序只要求路径存在且能载入）。
- 文件编码统一 UTF-8；CSV 允许带 BOM（程序以 `utf-8-sig` 读取）。

## 六、规模参考

| 谱库 | 条目数 | JSON 体积 | 载入后内存 |
| --- | --- | --- | --- |
| EDB 实验库 | 约 103 万 | 2.1 GB | — |
| ISDB 预测库 | 约 47 万 | 586 MB | — |
| 两者合计 | — | 约 2.7 GB | 峰值 5 GB 以上 |

两个 JSON 通过 `ujson` 一次性读入内存，因此建议机器内存 16 GB 以上。`examples/demo_db/` 里的精简库仅 96 KB，用于跑通流程。
