# 算法清单

引擎对每个候选分别计算以下 46 种谱相似度指标，每种指标的结果写入 `score__<name>` 列，同一指标下匹配上的碎片数写入 `npeak__<name>` 列。

> 本清单由输出表 `tournament_summary.csv` 的实际列名导出，顺序与列顺序一致。

| # | 指标标识 | 说明 |
|---|---|---|
| 1 | `dot_product_reverse` | 反向点积 |
| 2 | `mean_character` | 均值特征 |
| 3 | `squared_euclidean` | 平方欧氏距离 |
| 4 | `ms_for_id_v1` | MSforID 方法（v1） |
| 5 | `symmetric_chi_squared` | 对称卡方 |
| 6 | `weighted_dot_product` | 加权点积 |
| 7 | `chebyshev` | 切比雪夫距离 |
| 8 | `pearson_correlation` | 皮尔逊相关 |
| 9 | `penrose_shape` | Penrose 形状相似度 |
| 10 | `euclidean` | 欧氏距离 |
| 11 | `peak_percentage` | 峰占比 |
| 12 | `bhattacharya_1` | Bhattacharya 距离（1） |
| 13 | `bhattacharya_2` | Bhattacharya 距离（2） |
| 14 | `baroni_urbani_buser` | Baroni-Urbani-Buser 系数 |
| 15 | `squared_chord` | 平方弦距离 |
| 16 | `fidelity` | 保真度 |
| 17 | `dot_product` | 点积 |
| 18 | `modified_cosine` | 修正余弦 |
| 19 | `neutral_loss` | 中性丢失匹配 |
| 20 | `spectral_contrast_angle` | 谱对比角 |
| 21 | `unweighted_entropy` | 非加权熵相似度 |
| 22 | `harmonic_mean` | 调和平均 |
| 23 | `probabilistic_symmetric_chi_squared` | 概率对称卡方 |
| 24 | `vicis_symmetric_chi_squared_3` | VICIS 对称卡方（3） |
| 25 | `entropy` | 熵相似度 |
| 26 | `absolute_value` | 绝对值距离 |
| 27 | `intersection` | 交集 |
| 28 | `manhattan` | 曼哈顿距离 |
| 29 | `motyka` | Motyka 相似度 |
| 30 | `dice` | Dice 系数 |
| 31 | `lorentzian` | Lorentzian 距离 |
| 32 | `clark` | Clark 距离 |
| 33 | `improved_similarity` | 改进相似度 |
| 34 | `matusita` | Matusita 距离 |
| 35 | `ruzicka` | Ruzicka 相似度 |
| 36 | `jaccard` | Jaccard 系数 |
| 37 | `roberts` | Roberts 距离 |
| 38 | `penrose_size` | Penrose 尺寸相似度 |
| 39 | `ms_for_id` | MSforID |
| 40 | `hellinger` | Hellinger 距离 |
| 41 | `avg_l` | 平均 L 距离 |
| 42 | `whittaker_index_of_association` | Whittaker 关联指数 |
| 43 | `inner_product` | 内积 |
| 44 | `canberra` | Canberra 距离 |
| 45 | `wave_hedges` | Wave-Hedges 距离 |
| 46 | `divergence` | 散度 |

## 说明

- 表中「说明」列为各指标的通用名称释义，具体实现细节请参阅源码。
- 部分指标实现参考了 `matchms`（如 `modified_cosine`、`neutral_loss`）与 `ms_entropy`（FlashEntropySearch 相关熵类指标），向量距离类指标参考了 SciPy 的距离函数。
- 本项目的原创部分在于：将 46 路指标同时用于候选打分，再以锦标赛机制（协同得分提名 + 三轨决选）收敛到最终结果。
