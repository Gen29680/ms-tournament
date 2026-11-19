#!/usr/bin/env python3
import sys
import os
import multiprocessing
import argparse
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
from matchms import Spectrum
import gc  

# 移除暴力静默，恢复系统健康报警能力
_original_argv = sys.argv.copy()
sys.argv = [sys.argv[0]] 

try:
    from my_packages import functions_new, ms2tools_new
except ImportError:
    print("Error: 找不到 my_packages 模块。")
    sys.exit(1)
finally:
    sys.argv = _original_argv 

# 全局变量
GLOBAL_E2_JS = None
GLOBAL_I2_JS = None

METHODS_46 = [
    "dot_product_reverse","mean_character","squared_euclidean","ms_for_id_v1",
    "symmetric_chi_squared","weighted_dot_product","chebyshev","pearson_correlation",
    "penrose_shape","euclidean","peak_percentage","bhattacharya_1","bhattacharya_2",
    "baroni_urbani_buser","squared_chord","fidelity","dot_product","modified_cosine",
    "neutral_loss","spectral_contrast_angle","unweighted_entropy","harmonic_mean",
    "probabilistic_symmetric_chi_squared","vicis_symmetric_chi_squared_3","entropy",
    "absolute_value","intersection","manhattan","motyka","dice","lorentzian","clark",
    "improved_similarity","matusita","ruzicka","jaccard","roberts","penrose_size",
    "ms_for_id","hellinger","avg_l","whittaker_index_of_association","inner_product",
    "canberra","wave_hedges","divergence"
]

def parse_args():
    p = argparse.ArgumentParser(
        description="MS2 Tournament V6.8 - Ultimate Speed & Clean Edition",
        formatter_class=argparse.RawTextHelpFormatter # 允许 help 信息换行
    )
    
    # ------------------- 数据输入输出参数 (I/O) -------------------
    p.add_argument('--query_mgf', required=True,
                   help='[必填] 待鉴定的实验图谱文件(MGF格式)，包含采集到的二级质谱图。')
    p.add_argument('--edbms1_csv', required=True,
                   help='[必填] 实验标准数据库(EDB)的MS1索引文件(CSV格式)，用于根据前体离子进行初筛。')
    p.add_argument('--isms1_csv', required=True,
                   help='[必填] 内部/模拟数据库(ISDB)的MS1索引文件(CSV格式)，用于根据前体离子进行初筛。')
    p.add_argument('--edbms2_json', required=True,
                   help='[必填] 实验标准数据库(EDB)的MS2碎片图谱库(JSON格式)，用于进入二级深度比对。')
    p.add_argument('--isms2_json', required=True,
                   help='[必填] 内部/模拟数据库(ISDB)的MS2碎片图谱库(JSON格式)，用于进入二级深度比对。')
    p.add_argument('--output', default='ms_tournament_results',
                   help='[可选] 输出结果的根目录名称。默认为 "ms_tournament_results"。\n'
                        '脚本会在此目录下为每个成功匹配的 Query 动态创建专属子文件夹。')

    # ------------------- 质谱匹配阈值参数 (Thresholds) -------------------
    p.add_argument('--ppm_tol', type=float, default=10.0,
                   help='[可选] MS1前体离子质量容差 (ppm)，默认 10.0。\n'
                        '只有数据库中目标化合物的质量落在 Query m/z ± 该容差范围内，才会进入二级比对。')
    p.add_argument('--ms2_ppm_strict', type=float, default=20.0,
                   help='[可选] MS2碎片匹配严格容差 (ppm)，默认 20.0。\n'
                        '用于判定查询图谱和参考图谱的碎片峰是否为同一个峰。')
    p.add_argument('--noise_threshold', type=float, default=0.01,
                   help='[可选] 相对强度噪声剔除阈值，默认 0.01 (即 1%%)。\n'
                        '自动剔除相对于基峰强度低于此阈值的噪声峰，可大幅提升计算速度。')
    p.add_argument('--min_peaks', type=int, default=3,
                   help='[可选] 最少匹配碎片峰数量，默认 3。\n'
                        'Query和Reference之间匹配上的碎片峰如果少于该值，将被直接淘汰。')

    # ------------------- 系统控制与策略参数 (System & Strategy) -------------------
    p.add_argument('--cpus', type=int, default=8,
                   help='[可选] 并行处理的CPU核心/线程数，默认 8。\n'
                        '数值越高理论并发越高，但需根据物理机内存容量合理配置。')
    p.add_argument('--compare_all_is_energies', action='store_true',
                   help='[可选标志位] ISDB全碰撞能量比对开关。\n'
                        '触发时，将测试模拟库中所有碰撞能量(e0,e1,e2)并自动提取最高分；\n'
                        '不触发时，默认只比对第一种能量以节省时间。')
    p.add_argument('--no_resume', action='store_true',
                   help='[可选标志位] 禁用断点续传。\n'
                        '强制忽略已有文件，从头开始重新计算所有数据。')
    
    return p.parse_args()

def clean_spectrum(spec_arr, threshold):
    if spec_arr is None or len(spec_arr) == 0: return spec_arr
    spec_arr = np.array(spec_arr)
    max_i = np.max(spec_arr[:, 1])
    if max_i > 0:
        mask = spec_arr[:, 1] >= (max_i * threshold)
        return spec_arr[mask]
    return spec_arr

def make_mirror_plot(spec_a, spec_b, out_path, title=None):
    mz1, i1 = spec_a[:,0], spec_a[:,1]
    mz2, i2 = spec_b[:,0], spec_b[:,1]
    i1n = i1 / i1.max() if i1.max() > 0 else i1
    i2n = i2 / i2.max() if i2.max() > 0 else i2
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.vlines(mz1, 0, i1n, colors='tab:blue', linewidth=1.5, alpha=0.8, label='Query (Cleaned)')
    ax.vlines(mz2, 0, -i2n, colors='tab:orange', linewidth=1.5, alpha=0.8, label='Reference (Cleaned)')
    ax.axhline(0, color='black', linewidth=1)
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlabel('m/z')
    ax.set_ylabel('Relative Intensity')
    if title: ax.set_title(title, fontsize=9, pad=10)
    plt.legend(loc='upper right')
    plt.tight_layout()
    fig.savefig(out_path, dpi=100) # 调低 DPI 以加速图片生成
    plt.close(fig)

def compare_pair(SPEC1, SPECTRUM1, SPEC2, SPECTRUM2, methods):
    res = {}
    for m in methods:
        try:
            s, n = functions_new.clac_spec_sim(SPEC1, SPECTRUM1, SPEC2, SPECTRUM2, m)
        except:
            s, n = np.nan, np.nan
        res[f'score__{m}'] = s
        res[f'npeak__{m}'] = n
    return res

def run_one_match(match_entry, q_arr, pm, qid, global_methods, all_energy, noise_thresh, global_e2, global_i2):
    mid = str(match_entry['match_id'])
    mtype = match_entry['match_type']
    out = match_entry.copy()
    
    try:
        q_arr_clean = clean_spectrum(q_arr, noise_thresh)
        if len(q_arr_clean) < 3: return out
        q_obj = Spectrum(mz=q_arr_clean[:,0].astype(float), intensities=q_arr_clean[:,1].astype(float), metadata={"precursor_mz": pm, "id": qid})
        
        t_mz_val = match_entry.get('target_mz')
        local_methods = global_methods.copy()
        if pd.isna(t_mz_val) or t_mz_val is None:
            t_mz = pm  
            local_methods = [m for m in local_methods if m not in ["modified_cosine", "neutral_loss"]]
        else:
            t_mz = float(t_mz_val)

        if mtype == 'EDB':
            s2, obj2 = functions_new.get_spectra_from_einfo(global_e2, mid)
            if s2 is not None:
                s2_clean = clean_spectrum(s2, noise_thresh)
                if len(s2_clean) >= 3:
                    obj2_clean = Spectrum(mz=s2_clean[:,0].astype(float), intensities=s2_clean[:,1].astype(float), metadata={"precursor_mz": t_mz})
                    res = compare_pair(q_arr_clean, q_obj, s2_clean, obj2_clean, local_methods)
                    out.update(res); out['energy'] = 'exp'
        else:
            e_data = functions_new.get_spectra_from_isinfo(global_i2, mid)
            cands = [('e0', e_data[0], e_data[1]), ('e1', e_data[2], e_data[3]), ('e2', e_data[4], e_data[5])]
            best_s, best_res = -1, None
            for tag, s2, obj2 in (cands if all_energy else [cands[0]]):
                if s2 is None: continue
                s2_clean = clean_spectrum(s2, noise_thresh)
                if len(s2_clean) < 3: continue
                obj2_clean = Spectrum(mz=s2_clean[:,0].astype(float), intensities=s2_clean[:,1].astype(float), metadata={"precursor_mz": t_mz})
                res = compare_pair(q_arr_clean, q_obj, s2_clean, obj2_clean, local_methods)
                curr_max = np.nanmax([v for k,v in res.items() if 'score__' in k] or [0])
                if curr_max > best_s:
                    best_s = curr_max; best_res = res; best_res['_tag'] = tag
            if best_res: out.update(best_res); out['energy'] = best_res['_tag']
    except: pass
    return out

def pick_tournament_v6(df, min_peaks):
    if df.empty: return pd.DataFrame(), pd.DataFrame()
    nominations = []
    for m in METHODS_46:
        s_col, n_col = f'score__{m}' , f'npeak__{m}'
        if s_col in df.columns:
            df_sub = df[df[n_col] >= min_peaks].copy()
            if df_sub.empty: continue
            df_sub['synergy_metric'] = df_sub[s_col] * np.log1p(df_sub[n_col])
            top3 = df_sub.nlargest(3, 'synergy_metric').copy()
            for i, (_, row) in enumerate(top3.iterrows(), 1):
                nominations.append({
                    'nominating_method': m, 'rank_in_model': i,
                    'match_id': row['match_id'], 'match_type': row['match_type'],
                    'name': row.get('name'), 'formula': row.get('formula'),       
                    'smiles': row.get('smiles'), 'inchi': row.get('inchi'),           
                    'query_mz': row.get('query_mz'), 'target_mz': row.get('target_mz'),
                    'score': row[s_col], 'npeak': row[n_col],
                    'synergy_score': row['synergy_metric'], 'energy': row.get('energy', 'exp')
                })
    if not nominations: return pd.DataFrame(), pd.DataFrame()
    nom_df = pd.DataFrame(nominations)
    
    pool = df[df['match_id'].isin(nom_df['match_id'].unique())].copy()
    s_cols = [f'score__{m}' for m in METHODS_46]
    n_cols = [f'npeak__{m}' for m in METHODS_46]
    pool['global_max_score'] = pool[s_cols].max(axis=1)
    pool['global_max_peaks'] = pool[n_cols].max(axis=1)
    
    s_min, s_max = pool['global_max_score'].min(), pool['global_max_score'].max()
    p_min, p_max = pool['global_max_peaks'].min(), pool['global_max_peaks'].max()
    
    score_norm = (pool['global_max_score'] - s_min) / (s_max - s_min + 1e-6)
    peak_norm = (pool['global_max_peaks'] - p_min) / (p_max - p_min + 1e-6)
    pool['final_combined'] = (0.4 * score_norm) + (0.6 * peak_norm)
    
    r1 = pool.nlargest(10, 'global_max_score').drop_duplicates(subset=['match_id']).head(3).copy()
    r1['selection_reason'] = 'Final_Score'
    r2 = pool.nlargest(10, 'global_max_peaks').drop_duplicates(subset=['match_id']).head(3).copy()
    r2['selection_reason'] = 'Final_Peaks'
    r3 = pool.nlargest(10, 'final_combined').drop_duplicates(subset=['match_id']).head(3).copy()
    r3['selection_reason'] = 'Final_Combined'
    
    return pd.concat([r1, r2, r3]).drop_duplicates(subset=['match_id', 'selection_reason']), nom_df

# ----------------- 核心重构：将单个 Query 打包，外层并行 -----------------
def process_single_query(qid, qinfo_dict, e1_df, i1_df, args, global_e2, global_i2):
    pm = float(qinfo_dict['pepmass'])
    tol = pm * args.ppm_tol / 1e6
    matches = []
    
    # MS1 快速初筛
    for df, mtype in [(e1_df, 'EDB'), (i1_df, 'ISDB')]:
        mz_col = 'pepmass' if mtype == 'EDB' else 'm+h'
        if mz_col not in df.columns: continue
        hits = df[(df[mz_col] >= pm - tol) & (df[mz_col] <= pm + tol)]
        for _, r in hits.iterrows():
            if pd.notna(r['id']): 
                matches.append({
                    'match_id': r['id'], 'match_type': mtype, 'query_mz': pm,
                    'target_mz': r.get(mz_col), 'name': r.get('name', r.get('compound_name', 'N/A')), 
                    'formula': r.get('formula', 'N/A'), 'smiles': r.get('smiles', 'N/A'), 'inchi': r.get('inchi', 'N/A')
                })
    if not matches: return None
    
    q_spec = qinfo_dict['spectra']
    q_arr = np.array(q_spec)
    
    # 对所有 matches 计算打分
    results = [
        run_one_match(m, q_arr, pm, qid, METHODS_46, args.compare_all_is_energies, args.noise_threshold, global_e2, global_i2)
        for m in matches
    ]
    df_res = pd.DataFrame([r for r in results if 'match_id' in r])
    if df_res.empty: return None
    
    final_9, nom_df = pick_tournament_v6(df_res, args.min_peaks)
    
    # ======== 【核心修改】确保只有产生实质结果，才开始建立文件夹 ========
    if nom_df.empty or final_9.empty: 
        return None
    
    # 动态建立专属文件夹
    q_dir = os.path.join(args.output, f"Query_{qid}")
    os.makedirs(q_dir, exist_ok=True)
    
    legacy_syn_file = os.path.join(q_dir, f"Query_{qid}_nominations_synergy.csv")
    final_res_file = os.path.join(q_dir, f"Query_{qid}_final_results.csv")
    nom_df.to_csv(legacy_syn_file, index=False)
    
    final_9['matched_query_mz'] = ""
    final_9['matched_target_mz'] = ""
    final_9['extraction_track'] = "" 
    
    q_arr_clean = clean_spectrum(q_arr, args.noise_threshold)
    
    for idx, row in final_9.iterrows():
        try:
            mid, mt, reason = row['match_id'], row['match_type'], row['selection_reason']
            if mt == 'EDB': ref_s, _ = functions_new.get_spectra_from_einfo(global_e2, mid)
            else:
                ref_data = functions_new.get_spectra_from_isinfo(global_i2, mid)
                ref_s = ref_data[{'e0':0, 'e1':2, 'e2':4}.get(row['energy'], 0)]
            
            r_arr_clean = clean_spectrum(np.array(ref_s), args.noise_threshold)
            matched_q, matched_t = [], []
            track_status = "Strict_20ppm"
            
            for q_mz in q_arr_clean[:, 0]:
                dynamic_tol_strict = max(q_mz * args.ms2_ppm_strict / 1e6, 0.01)
                diffs = np.abs(r_arr_clean[:, 0] - q_mz)
                min_idx = np.argmin(diffs)
                if diffs[min_idx] <= dynamic_tol_strict:
                    matched_q.append(f"{q_mz:.4f}")
                    matched_t.append(f"{r_arr_clean[min_idx, 0]:.4f}")
            
            if not matched_q and row['global_max_peaks'] > 0:
                track_status = "Fallback_100ppm"
                for q_mz in q_arr_clean[:, 0]:
                    dynamic_tol_fallback = max(q_mz * (args.ms2_ppm_strict * 5) / 1e6, 0.02)
                    diffs = np.abs(r_arr_clean[:, 0] - q_mz)
                    min_idx = np.argmin(diffs)
                    if diffs[min_idx] <= dynamic_tol_fallback:
                        matched_q.append(f"{q_mz:.4f}")
                        matched_t.append(f"{r_arr_clean[min_idx, 0]:.4f}")
            
            final_9.at[idx, 'matched_query_mz'] = "; ".join(matched_q)
            final_9.at[idx, 'matched_target_mz'] = "; ".join(matched_t)
            final_9.at[idx, 'extraction_track'] = track_status

            title = f"[{reason}] ID:{mid} | Score:{row['global_max_score']:.3f} | Peaks:{int(row['global_max_peaks'])}"
            make_mirror_plot(q_arr_clean, r_arr_clean, os.path.join(q_dir, f"{reason}_{mid}.png"), title)
        except:
            final_9.at[idx, 'extraction_track'] = "ERROR"
            
    final_9['query_id'] = qid
    final_9.to_csv(final_res_file, index=False)
    
    return final_9


def main():
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)
    
    print(f"正在加载并预处理数据库... (CPU: {args.cpus} 核)")
    qinfo = functions_new.spectra_process(args.query_mgf)
    e1_df = pd.read_csv(args.edbms1_csv, low_memory=False)
    i1_df = pd.read_csv(args.isms1_csv, low_memory=False)
    
    # 提前转换数据类型，避免多进程内耗
    for df, col in [(e1_df, 'pepmass'), (i1_df, 'm+h')]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
            
    GLOBAL_E2_JS = functions_new.json_load(args.edbms2_json)
    GLOBAL_I2_JS = functions_new.json_load(args.isms2_json)

    # 预解构 Query 信息
    qinfo_tasks = {}
    for qid in qinfo.keys():
        q_spec, _ = functions_new.get_spectra_from_info(qinfo, qid)
        qinfo_tasks[qid] = {'pepmass': qinfo[qid]['pepmass'], 'spectra': q_spec}

    print(f"开始并行计算 {len(qinfo_tasks)} 个 Query 的锦标赛比对...")
    
    # 使用后端多线程或多进程进行大粒度并行
    global_summary = Parallel(n_jobs=args.cpus, backend='threading')(
        delayed(process_single_query)(qid, qinfo_tasks[qid], e1_df, i1_df, args, GLOBAL_E2_JS, GLOBAL_I2_JS)
        for qid in qinfo_tasks.keys()
    )
    
    # 过滤掉为 None 的空结果（即那些被我们淘汰，没有生成文件夹的 Query）
    global_summary = [res for res in global_summary if res is not None and not res.empty]

    if global_summary:
        cols = list(global_summary[0].columns)
        if 'query_id' in cols:
            cols.insert(0, cols.pop(cols.index('query_id')))
        pd.concat(global_summary)[cols].to_csv(os.path.join(args.output, 'tournament_summary.csv'), index=False)
        print(f"🎉 全部运行完毕！共成功鉴定 {len(global_summary)} 个 Query，无无效空文件夹产生。")
    else:
        print("⚠️ 运行结束，未找到任何符合条件的匹配结果。")

if __name__ == '__main__':
    if sys.platform != 'win32':
        try: multiprocessing.set_start_method('fork', force=True)
        except RuntimeError: pass
    main()