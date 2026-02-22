import os
import ast
import re
import numpy as np
import pandas as pd
from typing import List, Tuple, Dict, Optional
from .parsers.utils import read_table_auto
from sklearn.isotonic import spearmanr
from sklearn.metrics import mean_squared_error
from sklearn.metrics import (
    mean_squared_error, roc_auc_score, precision_recall_curve, auc,
    accuracy_score, precision_score, recall_score, f1_score,
    confusion_matrix, cohen_kappa_score, brier_score_loss
)
from .realbench.utils import process_variant_data
import scipy.io
from scipy.sparse import coo_matrix, csr_matrix
import sys
import subprocess



NA_STRINGS = {"", "NA", "NaN", "nan", "NA|NA", "nan|nan", "<NA>", "<NA>|<NA>"}
_region_re = re.compile(r"^chr[^:]+:(\d+)-(\d+)$", re.IGNORECASE)

def _infer_binsize_mode_from_cna(path: str, nrows: int = 5000) -> Optional[int]:
    df = pd.read_csv(path, sep=None, engine="python", nrows=nrows)
    region = df["region"] if "region" in df.columns else df.iloc[:, 0]

    lens = []
    for x in region.dropna().astype(str):
        m = _region_re.match(x.strip())
        if not m:
            continue
        start, end = int(m.group(1)), int(m.group(2))
        if end > start:
            lens.append(end - start+1)

    if not lens:
        return None

    vc = pd.Series(lens).value_counts()

    top = vc[vc == vc.iloc[0]].index.min()
    return int(top)

def check_binsize(tool_cna_files: List[str], tool_names: List[str], strict: bool = False) -> Dict[str, Optional[int]]:
    bs = {name: _infer_binsize_mode_from_cna(path) for path, name in zip(tool_cna_files, tool_names)}
    uniq = sorted({v for v in bs.values() if v is not None})

    if len(uniq) <= 1:
        print(f"[bin size check] pass: {uniq[0]}" if uniq else "ERROR: can not get bin size")
        return bs

    msg = "[bin size check] not pass: \n" + "\n".join([f"  {k}: {v}" for k, v in bs.items()]) + f"\n  unique={uniq}"
    if strict:
        raise ValueError(msg)
    print(msg)
    return bs

def align(df_pred: pd.DataFrame, df_true: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    common_idx = df_pred.index.intersection(df_true.index)
    common_col = df_pred.columns.intersection(df_true.columns)
    return df_pred.loc[common_idx, common_col], df_true.loc[common_idx, common_col]


def flatten_pair(a: pd.DataFrame, b: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    xa, xb = a.values.flatten().astype(float), b.values.flatten().astype(float)
    mask = ~np.isnan(xa) & ~np.isnan(xb)
    return xa[mask], xb[mask]

def to_numeric(df: pd.DataFrame) -> pd.DataFrame:
    return df.apply(pd.to_numeric, errors="coerce")


def safe_rmse(a: pd.DataFrame, b: pd.DataFrame) -> float:
    xa, xb = flatten_pair(a, b)
    if xa.size == 0:
        return float("nan")
    return float(np.sqrt(mean_squared_error(xb, xa)))


def safe_spearman(a: pd.DataFrame, b: pd.DataFrame) -> float:
    xa, xb =flatten_pair(a, b)
    if xa.size == 0:
        return float("nan")
    r, _ = spearmanr(xb, xa)
    return float(r)


def safe_acc(a: pd.DataFrame, b: pd.DataFrame) -> float:
    va, vb = a.values, b.values
    mask = ~pd.isna(va) & ~pd.isna(vb)
    n = mask.sum()
    if n == 0:
        return float("nan")
    return float((va[mask] == vb[mask]).mean())


def to_numeric(df: pd.DataFrame) -> pd.DataFrame:
    return df.apply(pd.to_numeric, errors="coerce")


def safe_acc_combined(pred_df: pd.DataFrame, truth_df: pd.DataFrame) -> float:
    a_str = pred_df.astype("string")
    b_str = truth_df.astype("string")

    mask = ~(a_str.isna() | b_str.isna())

    valid = mask.to_numpy()

    a_vals = a_str.to_numpy()
    b_vals = b_str.to_numpy()

    a_valid = a_vals[valid]
    b_valid = b_vals[valid]

    if len(a_valid) == 0:
        return float("nan")

    
    return float((a_valid == b_valid).mean())

def split_from_combined(df: pd.DataFrame) -> pd.DataFrame:

    df = df.astype(str)

    s0 = df.apply(lambda x: x.str.split("|").str[0]).apply(pd.to_numeric, errors="coerce")
    s1 = df.apply(lambda x: x.str.split("|").str[1]).apply(pd.to_numeric, errors="coerce")
    return s0 , s1



def evaluate_haplotype_predictions(
    pred_df: pd.DataFrame, truth_df: pd.DataFrame, haplotype: str
) -> Tuple[float, float, float]:

    pred_df, truth_df = align(pred_df, truth_df)

    if haplotype == "combined":

        pred_a, pred_b = split_from_combined(pred_df)
        truth_a,truth_b = split_from_combined(truth_df)
        print("complete ")
        rmse = (safe_rmse(pred_a, truth_a) + safe_rmse(pred_b, truth_b)) /2
        scc  = ( safe_spearman(pred_a, truth_a) + safe_spearman(pred_b, truth_b))/2

        acc  = safe_acc_combined(pred_df, truth_df)
    else:

        pred_num = to_numeric(pred_df)
        truth_num = to_numeric(truth_df)
        rmse = safe_rmse(pred_num, truth_num)
        scc  = safe_spearman(pred_num, truth_num)
        acc  = safe_acc(pred_num, truth_num)

    return rmse, scc, acc

def calc_bin_metrics(gdf: pd.DataFrame, pdf: pd.DataFrame) -> Dict[str, Optional[float]]:
        y_true = gdf.values.flatten()
        y_pred = pdf.values.flatten()
        metrics = {}
        metrics["Accuracy"]  = accuracy_score(y_true, y_pred)
        metrics["F1"]        = f1_score(y_true, y_pred, zero_division=0)
        metrics["Brier"]     = brier_score_loss(y_true, y_pred)

        if len(set(y_true)) > 1 and len(set(y_pred)) > 1:
            labels = sorted(set(y_true).union(set(y_pred)))
            precision, recall, _ = precision_recall_curve(y_true, y_pred)
            metrics["AUROC"] = roc_auc_score(y_true, y_pred)
            metrics["AUPRC"] = auc(recall, precision)
            tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=labels).ravel()
            metrics["Sensitivity"] = tp / (tp + fn) if (tp + fn) else None
            metrics["Specificity"] = tn / (tn + fp) if (tn + fp) else None
            metrics["PPV"] = tp / (tp + fp) if (tp + fp) else None
            metrics["NPV"] = tn / (tn + fn) if (tn + fn) else None

        else:
            metrics["AUROC"] = None
            metrics["AUPRC"] = None
            metrics["Kappa"]  = None
            metrics["Sensitivity"] = None
            metrics["Specificity"] = None
            metrics["PPV"] = None
            metrics["NPV"] = None
            

        return metrics 

def cat_fn_from_cond(condition: str):
    if condition.startswith(">="):
        thr = float(condition[2:])
        return lambda v: int(float(v) >= thr)
    if condition.startswith("="):
        thr = float(condition[1:])
        return lambda v: int(float(v) == thr)
    raise ValueError(f"no support condition: {condition}")

def metrics_per_clone(gdf: pd.DataFrame, pdf: pd.DataFrame) -> Dict[str, Dict]:
    out = {}

    clones = pdf.columns[1:].str.split("_").str[0].unique()
    for clone in clones:
        gsub = gdf.loc[:, gdf.columns.str.startswith(clone)]
        psub = pdf.loc[:, pdf.columns.str.startswith(clone)]
        out[clone] = calc_bin_metrics(gsub, psub)
    return out

def categorize_and_save(df_truth: pd.DataFrame, df_pred: pd.DataFrame,
                             out_dir: str, tool: str, condition: str, haplo_type: str):
        cat = cat_fn_from_cond(condition)
        gt_c = df_truth.copy()
        pd_c = df_pred.copy()
        for col in df_truth.columns[1:]:
            gt_c[col] = df_truth[col].apply(cat)
            pd_c[col] = df_pred[col].apply(cat)
        gt_c.to_csv(os.path.join(out_dir, f"ground_truth_{haplo_type}.csv"), index=False)
        pd_c.to_csv(os.path.join(out_dir, f"{tool}_predict_{haplo_type}.csv"), index=False)

def process_folder_for_metrics_clone(folder: str, tool: str, hap_list: List[str]) -> Dict[str, Dict[str, Dict]]:
    results = {}
    for h in hap_list:
        gpath = os.path.join(folder, f"ground_truth_{h}.csv")
        ppath = os.path.join(folder, f"{tool}_predict_{h}.csv")
        if not (os.path.exists(gpath) and os.path.exists(ppath)):
            continue
        gdf, pdf = pd.read_csv(gpath), pd.read_csv(ppath)
        gdf, pdf = align(gdf, pdf)
        results[h] = metrics_per_clone(gdf, pdf)
    return results

def process_folder_for_metrics(folder: str, tool: str, hap_list: List[str]) -> Dict[str, Dict[str, Dict]]:
    results = {}
    for h in hap_list:
        gpath = os.path.join(folder, f"ground_truth_{h}.csv")
        ppath = os.path.join(folder, f"{tool}_predict_{h}.csv")
        if not (os.path.exists(gpath) and os.path.exists(ppath)):
            continue
        gdf, pdf = pd.read_csv(gpath,index_col=0), pd.read_csv(ppath,index_col=0)
        gdf, pdf = align(gdf, pdf)
        results[h] = calc_bin_metrics(gdf, pdf)
    return results


def split_region(region_str):
    try:
        chrom, coords = region_str.split(":")
        start, end = map(int, coords.split("-"))
        return chrom, start, end
    except Exception as e:
        raise ValueError(f"Invalid segment format: {region_str}") from e


def add_coords(df):
    tmp = df["region"].str.split(":", expand = True)
    chrom = tmp[0]
    coords = tmp[1].str.split("-", expand = True)
    
    coord_df = pd.DataFrame({
        "chrom": chrom,
        "start": coords[0].astype(int),
        "end":   coords[1].astype(int),
    }, index=df.index)

    out = pd.concat([df, coord_df], axis=1).copy()
    return out

def build_master_bins(df1, df2):
    bins = []

    for chrom in sorted(set(df1.chrom) | set(df2.chrom)):
        a = df1[df1.chrom == chrom]
        b = df2[df2.chrom == chrom]

        a_starts = set(a.start)
        a_ends   = set(a.end + 1)   
        
        b_starts = set(b.start)
        b_ends   = set(b.end + 1)

        boundaries = sorted(a_starts | a_ends | b_starts | b_ends)

        for i in range(len(boundaries) - 1):
            s = boundaries[i]
            e = boundaries[i+1]
            bins.append([chrom, s, e - 1]) 

    return pd.DataFrame(bins, columns=["chrom", "start", "end"])

def map_to_master(master, old_df):
    old_df = old_df.copy()
    value_cols = [c for c in old_df.columns if c not in ["region", "chrom", "start", "end"]]

    out = pd.concat(
        [
            master[["chrom", "start", "end", "region"]].copy(),
            pd.DataFrame(np.nan, index=master.index, columns=value_cols),
        ],
        axis=1,
    )

    for chrom in master["chrom"].unique():
        m = master[master.chrom == chrom]
        o = old_df[old_df.chrom == chrom]

        if o.empty or m.empty:
            continue

        o = o.sort_values("start")
        starts = o["start"].to_numpy()
        ends   = o["end"].to_numpy()

        m_starts = m["start"].to_numpy()
        m_ends   = m["end"].to_numpy()

        idxs = np.searchsorted(starts, m_starts, side="right") - 1

        valid = (idxs >= 0) & (idxs < len(starts))

        valid &= ends[idxs] >= m_ends

        master_idx = m.index[valid]
        old_idx    = o.index[idxs[valid]]

        out.loc[master_idx, value_cols] = old_df.loc[old_idx, value_cols].to_numpy()

    return out


def align_cna_bins(cna_df1,cna_df2):

    cna_df1 = add_coords(cna_df1)
    cna_df2 = add_coords(cna_df2)

    df1 = cna_df1[['chrom', 'start', 'end']]
    df2 = cna_df2[['chrom', 'start', 'end']] 

    cna_df1 = cna_df1.sort_values(["chrom", "start"]).reset_index(drop=True)
    cna_df2 = cna_df2.sort_values(["chrom", "start"]).reset_index(drop=True)

    master = build_master_bins(df1, df2)
    master["region"] = master.apply(
        lambda r: f"{r.chrom}:{r.start}-{r.end}", axis=1
    )

    cna_df1_new = map_to_master(master,cna_df1)
    cna_df1_new.drop(columns=['chrom', 'start', 'end'],inplace=True)

    cna_df2_new = map_to_master(master,cna_df2)
    cna_df2_new.drop(columns=['chrom', 'start', 'end'],inplace=True)


    return cna_df1_new, cna_df2_new

def fast_mode_num(arr):
    arr = arr[~np.isnan(arr)]
    if arr.size == 0:
        return None
    vals, counts = np.unique(arr, return_counts=True)
    return vals[np.argmax(counts)]   

def fast_mode_any(x):
    x = np.asarray(x, dtype=object)

    mask = ~pd.isna(x)
    x = x[mask]
    if x.size == 0:
        return None

    xs = np.array([str(v).strip() for v in x], dtype=object)
    xs = xs[~np.isin(xs, list(NA_STRINGS))]
    if xs.size == 0:
        return None

    vals, cnt = np.unique(xs, return_counts=True)
    return vals[cnt.argmax()]


def read_and_drop_empty(path):
    df = pd.read_csv(path)

    df = df.dropna(axis=0, how='all')
    df = df.dropna(axis=1, how='all')

    return df

def evaluate_NA_ratio(df: pd.DataFrame, region_col: str = "region"):

    cols = [c for c in df.columns if c != region_col]
    if not cols:
        return 0.0, 0, 0

    sub = df[cols]

    na_mask = sub.isna()

    obj_cols = sub.select_dtypes(include=["object", "string"]).columns
    if len(obj_cols) > 0:
        empty_mask = sub[obj_cols].apply(lambda s: s.astype("string").str.strip().eq(""))
        na_mask.loc[:, obj_cols] = na_mask.loc[:, obj_cols] | empty_mask

    na_count = int(na_mask.to_numpy().sum())
    total = int(sub.shape[0] * sub.shape[1])
    na_ratio = na_count / total if total else 0.0
    return na_ratio, na_count, total

def _get_ref_path(filename):
        return os.path.join((os.path.dirname(__file__)), "ref", filename)

def load_centromeres(genome="hg38"):
        print(f"[hcbench] Loading {genome} centromere data")
        df = pd.read_csv(_get_ref_path(f"{genome}_centromeres.csv"))
        # print(f"[hcbench] Loading centromere data from {data_path}...")
        return df

def _classify_size(start, end, small_size, mid_size):
        size = int(end) - int(start) + 1
        if size < small_size:
            return "focal"
        elif size <= mid_size:
            return "medium"
        return "broad"

def _classify_cnv(val: str) -> str:
        try:
            left, right = map(int, str(val).split("|"))
        except Exception:
            return "UNKNOWN"

        s = left + right
        if s == 0:
            return "DEL"
        elif s == 1:
            return "CNL_LOH"
        elif s == 2 and (left == 0 or right == 0):
            return "CNN_LOH"
        elif s > 2 and (left == 0 or right == 0):
            return "CNG_LOH"
        elif s > 2 and (left != 0 and right != 0):
            return "DUP"
        return "UNKNOWN"

def _classify_whole_arm_events(
        val: str, chrom: str, start: int, end: int, centromeres_df: pd.DataFrame
    ) -> str:
        try:
            left, right = map(int, str(val).split("|"))
        except Exception:
            return "UNKNOWN"
        
        start = int(start)
        end = int(end)

        ch = str(chrom)
        if not ch.startswith("chr"):
            ch = "chr" + ch

        row = centromeres_df.loc[centromeres_df["chrom"] == ch]
        if row.empty:
            return _classify_cnv(val)
        row = row.iloc[0]

        p_end = int(row["p_arm_end"])
        q_start = int(row["q_arm_start"])
        chr_len = int(row["chrom_length"])

        covers_p = (start <= 1) and (end >= p_end)
        covers_q = (start <= q_start) and (end >= chr_len)

        # if (left == 0 or right == 0) and ((start <= 1) and (end >= chr_len)):
        #     return "WCL"  # Whole-arm loss
        # if (left > 1 or right > 1) and ((start <= 1) and (end >= chr_len)) and (left != 0 and right != 0):
        #     return "WCD"  # Whole-arm gain
        
        if (left == 0 or right == 0) and (covers_p or covers_q):
            return "ARM-DEL"  # arm loss
        if (left > 1 or right > 1) and (covers_p or covers_q) and (left != 0 and right != 0):
            # print("ARM-DUP DEBUG:", ch, start, end, "p_end", p_end, "q_start", q_start, "chr_len", chr_len, "covers_p", covers_p, "covers_q", covers_q, "val", val)
            return "ARM-DUP"  # arm gain

        return _classify_cnv(val)

def _mark_cell_wgd(df, centromeres_df, THRESH=0.8, WGD_SET = {"DUP", "ARM-DUP", "WCD", "WGD"}) -> pd.DataFrame:
    def _mark(cell_df):
        len_sum = centromeres_df['chrom_length'].astype(int).sum()
        cell_df['len'] = cell_df['End'].astype(int) - cell_df['Start'].astype(int) +1

        mask = cell_df["type"].isin(WGD_SET)
        ratio = cell_df.loc[mask, "len"].sum() / len_sum
        # ratio = cell_df["type"].isin(WGD_SET).mean()
        if ratio >= THRESH:
            cell_df.loc[mask, "type"] = "WGD"
        return cell_df

    return df.groupby("cell", group_keys=False).apply(_mark)

def _mark_chrom_wcd_wcl(df,centromeres_df, THRESH=0.8, WCD_SET = {"DUP", "ARM-DUP", "WCD"}, WCL_SET = {"DEL", "ARM-DEL","WCL"}) -> pd.DataFrame:
    def _mark(g):
        g['len'] = g['End'].astype(int) - g['Start'].astype(int) +1

        chrom = g["Segment"].iloc[0].split(":")[0]
        ch = chrom if chrom.startswith("chr") else "chr" + chrom
        row = centromeres_df.loc[centromeres_df["chrom"] == ch]
        chr_len = int(row["chrom_length"].iloc[0])

        if (g["type"] == "WGD").all():
            return g

        mask_WCD = g["type"].isin(WCD_SET)
        ratio = g.loc[mask_WCD, "len"].sum() / chr_len
        if ratio >= THRESH:
            g.loc[mask_WCD, "type"] = "WCD"

        mask_WCL = g["type"].isin(WCL_SET)
        ratio = g.loc[mask_WCL, "len"].sum() / chr_len
        if ratio >= THRESH:
            g.loc[mask_WCL, "type"] = "WCL"

        return g

    return df.groupby(["cell", "Chrom"], group_keys=False).apply(_mark)


def _mark_arm_events(df, centromeres_df, THRESH=0.8, ARM_DUP_SET = {"DUP", "ARM-DUP"}, ARM_DEL_SET = {"DEL", "ARM-DEL"}) -> pd.DataFrame:
    def _mark(g):
        if g["type"].isin({"WGD", "WCD","WCL"}).all():
            return g
        # print("Processing group: ", g.head())
        chrom = g["Segment"].iloc[0].split(":")[0]
        ch = chrom if chrom.startswith("chr") else "chr" + chrom
        g["Start"] = g["Start"].astype(int)
        g["End"]   = g["End"].astype(int)

        row = centromeres_df.loc[centromeres_df["chrom"] == ch]
        if row.empty:
            return g

        row = row.iloc[0]
        p_end = int(row["p_arm_end"])
        q_start = int(row["q_arm_start"])
        chr_len = int(row["chrom_length"])

        # ---- p arm ----
        p_bins = g[g["End"] <= p_end]
        if len(p_bins) > 0 :
            covers_p_start = p_bins["Start"].min() <= 1
            covers_p_right = (p_bins["End"].max() >= p_end)
            sum_len = (p_bins['End'] - p_bins['Start'] +1) .sum()
            if covers_p_start and covers_p_right and sum_len >= (p_end) * THRESH:
                dup_ratio = p_bins["type"].isin(ARM_DUP_SET).mean()
                del_ratio = p_bins["type"].isin(ARM_DEL_SET).mean()

                if dup_ratio >= THRESH:
                    g.loc[p_bins.index, "type"] = "ARM-DUP"
                elif del_ratio >= THRESH:
                    g.loc[p_bins.index, "type"] = "ARM-DEL"

        # ---- q arm ----
        q_bins = g[g["Start"] >= q_start]
        if len(q_bins) > 0:
            covers_q_end = q_bins["End"].max() >= chr_len
            covers_q_start = q_bins["Start"].min() <= q_start
            sum_len = (q_bins['End'] - q_bins['Start'] +1) .sum()
            if covers_q_end and covers_q_start and sum_len >= (chr_len - q_start) * THRESH:
                dup_ratio = q_bins["type"].isin(ARM_DUP_SET).mean()
                del_ratio = q_bins["type"].isin(ARM_DEL_SET).mean()

                if dup_ratio >= THRESH:
                    g.loc[q_bins.index, "type"] = "ARM-DUP"
                elif del_ratio >= THRESH:
                    g.loc[q_bins.index, "type"] = "ARM-DEL"

        return g

    return df.groupby(["cell", "Chrom"], group_keys=False).apply(_mark)


def annotate_segments(
        df: pd.DataFrame,
        small_size: int = 3_000_000,
        mid_size: int = 10_000_000,
        genome: str = "hg38",
        threshold: float = 0.8,
        detect_cnv_type = True
    ) -> pd.DataFrame:

        centromeres_df = load_centromeres(genome)
        merged_list = []
        cur_bin = []

        def flush_segment(cur_bin, chrom, start, end, cur_val):
            """
            Finalize one CNV segment (bin-level rows stored in cur_bin)
            """
            if not cur_bin:
                return
            cur_df = pd.DataFrame(cur_bin)
            cur_df["Segment"] = f"{chrom}:{start}-{end}"
            cur_df["size"] = _classify_size(start, end, small_size, mid_size)
            if detect_cnv_type:
                cur_df["type"] = _classify_whole_arm_events(
                    cur_val, chrom, start, end, centromeres_df
                )
            merged_list.append(cur_df)
        

        for cell, row in df.iterrows():
            chrom = start = end = cur_val = None
            cur_bin = []
            # print(f"Processing cell: {cell}")

            for key, val in row.items():
                chrom_new = key.split(":")[0]
                start_new = key.split("-")[0].split(":")[-1]
                end_new = key.split("-")[-1]

                # -------- neutral / invalid bin --------
                if (
                    val in ("1|1", "nan|nan")
                    or pd.isna(val)
                ):
                    flush_segment(cur_bin, chrom, start, end, cur_val)
                    chrom = start = end = cur_val = None
                    cur_bin = []
                    continue

                # -------- determine if new segment --------
                new_segment = (
                    chrom is None
                    or chrom_new != chrom
                    or val != cur_val
                )

                if new_segment:
                    flush_segment(cur_bin, chrom, start, end, cur_val)
                    # update current segment info
                    chrom = chrom_new
                    start = start_new
                    end = end_new
                    cur_val = val
                    cur_bin = []
                else:
                    # extend current segment
                    end = end_new

                # -------- append bin --------
                cur_bin.append({
                    "cell": cell,
                    "Chrom": chrom,
                    "Start": start_new,
                    "End": end_new,
                    "value": cur_val
                })


            # -------- finalize last segment of cell --------
            flush_segment(cur_bin, chrom, start, end, cur_val)

        final_df = pd.concat(merged_list, ignore_index=True)

        if not detect_cnv_type:
            return final_df

        if threshold >= 1:
            return final_df

        WGD_SET = {"ARM-DUP","DUP",}
        WCD_SET = {"ARM-DUP","DUP"}
        WCL_SET = {"ARM-DEL","CNN_LOH","CNL_LOH","CNG_LOH","DEL"}

        # final_df = _mark_arm_events(final_df, centromeres_df, THRESH=threshold, ARM_DUP_SET=ARM_DUP_SET, ARM_DEL_SET=ARM_DEL_SET)

        final_df_wcd_wcl = _mark_chrom_wcd_wcl(final_df,centromeres_df, THRESH=threshold, WCD_SET=WCD_SET, WCL_SET=WCL_SET)

        final_df_wcd_wcl = final_df_wcd_wcl[final_df_wcd_wcl['type'].isin(["WCL", "WCD"])]
        final_df_wcd_wcl['size'] = 'broad'
        final_df_wcd_wcl.drop(columns=['len'],inplace=True)


        final_df_wgd = _mark_cell_wgd(final_df,centromeres_df, THRESH=threshold, WGD_SET=WGD_SET)
        final_df_wgd = final_df_wgd[final_df_wgd ['type'] == "WGD"]
        final_df_wgd['size'] = 'broad'
        final_df_wgd.drop(columns=['len'],inplace=True)

        final_df = pd.concat([final_df_wcd_wcl, final_df_wgd, final_df], ignore_index=True)



        return final_df 

def get_cell_profile_size(cn_df):

    groups = cn_df.T.groupby(list(cn_df.index), sort=False)

    cluster_id = groups.ngroup() 

    cluster_size = cluster_id.map(cluster_id.value_counts())

    result = pd.DataFrame({
        "cell": cn_df.columns,
        # "cluster_id": cluster_id.values,
        "cluster_size": cluster_size.values
    })

    result.set_index("cell",inplace = True)

    return result 

def get_cluster_size(cluster_df: pd.DataFrame) -> pd.DataFrame:
    cluster_id = cluster_df["clone_id"]

    if cluster_id.isna().any():
        cluster_id = cluster_id.copy()
        na_cells = cluster_id.isna()
        cluster_id.loc[na_cells] = ["singleton_" + str(i) for i in range(na_cells.sum())]

    counts = cluster_id.value_counts()
    cluster_size = cluster_id.map(counts)

    out = pd.DataFrame({"cluster_id": cluster_id, "cluster_size": cluster_size}, index=cluster_df.index)
    out.index.name = "cell"
    return out

def get_segment_overlap_ratio(gt_annotated_df, pred_annotated_df):
    sizes = ['focal', 'medium', 'broad']
    type_list = sorted(set(gt_annotated_df["type"].unique()) | set(pred_annotated_df["type"].unique()))
    cols = ['cell', 'region']

    out = []
    for s in sizes:
        for t in type_list:
            gt_s = gt_annotated_df.loc[(gt_annotated_df['size'].eq(s) & gt_annotated_df['type'].eq(t)),  cols].drop_duplicates()
            pred_s = pred_annotated_df.loc[(pred_annotated_df['size'].eq(s) & pred_annotated_df['type'].eq(t)), cols].drop_duplicates()
            
            if gt_s.empty or pred_s.empty:
                out.append((s, t, np.nan, np.nan, np.nan))
                continue

            for col in ['cell', 'region']:
                gt_s[col] = gt_s[col].astype('category')
                pred_s[col] = pred_s[col].astype('category')

            inter_n = gt_s.merge(pred_s, on=cols, how='inner').shape[0]
            union_n = pd.concat([gt_s, pred_s], ignore_index=True).drop_duplicates().shape[0]

            out.append((s, t,union_n, inter_n, inter_n / union_n if union_n else 0.0))

    res = pd.DataFrame(out, columns=['Size', 'Type','Sum_n', 'overlap_n', 'overlap_ratio'])

    return res

def get_segment_metric_both(gt_annotated_df, pred_annotated_df):
    sizes = ['focal', 'medium', 'broad']
    type_list = set(gt_annotated_df['type'].unique()).intersection(set(pred_annotated_df['type'].unique()))
    cols = ['cell', 'region']

    out = []
    for s in sizes:
        for t in type_list:
            gt_s = gt_annotated_df.loc[(gt_annotated_df['size'].eq(s) & gt_annotated_df['type'].eq(t)),  cols + ["value"]].drop_duplicates()
            pred_s = pred_annotated_df.loc[(pred_annotated_df['size'].eq(s) & pred_annotated_df['type'].eq(t)),  cols + ["value"]].drop_duplicates()

            for col in ['cell', 'region']:
                gt_s[col] = gt_s[col].astype('category')
                pred_s[col] = pred_s[col].astype('category')

            inter_df = gt_s.merge(pred_s, on=cols, how='inner')

            if inter_df.empty:
                out.append((s, t, np.nan, np.nan, np.nan))
                continue

            inter_df['gt_cn_a'] = inter_df['value_x'].str.split('|').str[0].astype(int)
            inter_df['gt_cn_b'] = inter_df['value_x'].str.split('|').str[1].astype(int)

            inter_df['tool_cn_a'] = inter_df['value_y'].str.split('|').str[0].astype(int)
            inter_df['tool_cn_b'] = inter_df['value_y'].str.split('|').str[1].astype(int)

            acc = (inter_df['value_x'] == inter_df['value_y']).mean()
            rmse = np.sqrt(mean_squared_error(inter_df['gt_cn_a'], inter_df['tool_cn_a']) + mean_squared_error(inter_df['gt_cn_b'], inter_df['tool_cn_b'])) /2
            # mse = ((inter_df['gt_cn_a'] - inter_df['tool_cn_a']) ** 2 + (inter_df['gt_cn_b'] - inter_df['tool_cn_b']) ** 2).mean()
            scc = (spearmanr(inter_df['gt_cn_a'], inter_df['tool_cn_a'])[0] + spearmanr(inter_df['gt_cn_b'], inter_df['tool_cn_b'])[0]) / 2

            out.append((s,t, acc, scc, rmse))

    res = pd.DataFrame(out, columns=['Size','Type', 'ACC', 'SCC', 'RMSE'])
    return res

def get_segment_metric(gt_annotated_df, pred_df):
    sizes = ['focal', 'medium', 'broad']
    type_list = gt_annotated_df['type'].unique()
    cols = ['cell', 'region']

    out = []
    for s in sizes:
        for t in type_list:
            gt_s = gt_annotated_df.loc[(gt_annotated_df['size'].eq(s) & gt_annotated_df['type'].eq(t)),  cols + ["value"]].drop_duplicates()
            pred_s = pred_df.drop_duplicates()

            for col in ['cell', 'region']:
                gt_s[col] = gt_s[col].astype('category')
                pred_s[col] = pred_s[col].astype('category')

            inter_df = gt_s.merge(pred_s, on=cols, how='inner')

            if inter_df.empty:
                out.append((s, t, np.nan, np.nan, np.nan))
                continue

            inter_df['gt_cn_a'] = inter_df['value_x'].str.split('|').str[0].astype(int)
            inter_df['gt_cn_b'] = inter_df['value_x'].str.split('|').str[1].astype(int)

            inter_df['tool_cn_a'] = inter_df['value_y'].str.split('|').str[0].astype(int)
            inter_df['tool_cn_b'] = inter_df['value_y'].str.split('|').str[1].astype(int)

            acc = (inter_df['value_x'] == inter_df['value_y']).mean()
            rmse = np.sqrt(mean_squared_error(inter_df['gt_cn_a'], inter_df['tool_cn_a']) + mean_squared_error(inter_df['gt_cn_b'], inter_df['tool_cn_b'])) /2
            # mse = ((inter_df['gt_cn_a'] - inter_df['tool_cn_a']) ** 2 + (inter_df['gt_cn_b'] - inter_df['tool_cn_b']) ** 2).mean()
            scc = (spearmanr(inter_df['gt_cn_a'], inter_df['tool_cn_a'])[0] + spearmanr(inter_df['gt_cn_b'], inter_df['tool_cn_b'])[0]) / 2

            out.append((s, t, acc, scc, rmse))

    res = pd.DataFrame(out, columns=['Size', 'Type', 'ACC', 'SCC', 'RMSE'])
    return res

def get_seg_cna_event_num(df):

    bin_level_cna_event_counts = (
        df
        .groupby(["size", "type"])
        .size()
        .reset_index(name="bin_level_count")
    )

    cell_seg_uniq = (
        df
        .groupby(["size", "type", "cell"])["Segment"]
        .nunique(dropna=True)
        .reset_index(name="seg_uniq_per_cell")
    )

    # Step 2: 对每个 (size, type) 把各 cell 的 unique 数量相加
    seg_level_cna_event_counts = (
        cell_seg_uniq
        .groupby(["size", "type"])["seg_uniq_per_cell"]
        .sum()
        .reset_index(name="seg_level_count")
    )

    merged_cna_event_count = bin_level_cna_event_counts.merge(seg_level_cna_event_counts, on=["size", "type"], how="outer")

    return merged_cna_event_count



def find_mirrored_pairs(df: pd.DataFrame, clone_cols):
    
    long_df = (
        df[clone_cols]
        .stack(dropna=False)
        .reset_index()
        .rename(columns={"level_0": "region", "level_1": "Clone", 0: "CNA"})
    )


    ab = long_df["CNA"].astype(str).str.split("|", n=1, expand=True)
    long_df["a"] = pd.to_numeric(ab[0], errors="coerce")
    long_df["b"] = pd.to_numeric(ab[1], errors="coerce")

   
    long_df = long_df.dropna(subset=["a", "b"])
    long_df = long_df[long_df["a"] != long_df["b"]]

    long_df["a_i"] = long_df["a"].astype(int)
    long_df["b_i"] = long_df["b"].astype(int)
    long_df["key"] = long_df["a_i"].astype(str) + "|" + long_df["b_i"].astype(str)
    long_df["rev_key"] = long_df["b_i"].astype(str) + "|" + long_df["a_i"].astype(str)

    m = long_df.merge(
        long_df,
        left_on=["region", "key"],
        right_on=["region", "rev_key"],
        suffixes=("_1", "_2"),
        how="inner",
    )

    m = m[m["Clone_1"] < m["Clone_2"]]

    out = m.rename(
        columns={
            "Clone_1": "Clone1",
            "Clone_2": "Clone2",
            "CNA_1": "Clone1_CNA",
            "CNA_2": "Clone2_CNA",
        }
    )[["region", "Clone1", "Clone2", "Clone1_CNA", "Clone2_CNA"]]

    return out.reset_index(drop=True)

def find_mirrored_clones(df: pd.DataFrame) -> pd.DataFrame:
    print("[hcbench] Detecting mirrored clone CNAs...")
    clone_cols = df.columns
    mirrored_rows = []

    # for seg, row in df.iterrows():
    #     for i in range(len(clone_cols)):
    #         for j in range(i + 1, len(clone_cols)):
    #             c1, c2 = row[clone_cols[i]], row[clone_cols[j]]
    #             if pd.isna(c1) or pd.isna(c2):
    #                 continue
    #             try:
    #                 hap1 = tuple(map(int, c1.split("|")))
    #                 hap2 = tuple(map(int, c2.split("|")))
    #             except ValueError:
    #                 continue
    #             if hap1 == hap2[::-1] and hap1[0] != hap1[1]:
    #                 mirrored_rows.append(
    #                     {
    #                         "region": seg,
    #                         "Clone1": clone_cols[i],
    #                         "Clone2": clone_cols[j],
    #                         "Clone1_CNA": c1,
    #                         "Clone2_CNA": c2,
    #                     }
    #                 )

    # mirrored_df = pd.DataFrame(mirrored_rows)
    mirrored_df = find_mirrored_pairs(df, clone_cols)
    print(f"[hcbench] Found {len(mirrored_df)} mirrored clone events.")
    return mirrored_df



def get_seg_mirror_subclonal(mirrored_df,
        small_size: int = 3_000_000,
        mid_size: int = 10_000_000
    ):
    need_cols = {"Chromosome", "Start", "End"}
    if not need_cols.issubset(mirrored_df.columns):
        if "region" not in mirrored_df.columns:
            raise ValueError("Input df must contain either ['Chromosome','Start','End'] or 'region'.")

        chrom = mirrored_df["region"].astype(str).str.split(":", n=1).str[0]
        pos = mirrored_df["region"].astype(str).str.split(":", n=1).str[1]

        start = pos.str.split("-", n=1).str[0]
        end = pos.str.split("-", n=1).str[1]

        mirrored_df["Chromosome"] = chrom
        mirrored_df["Start"] = pd.to_numeric(start, errors="raise").astype(int)
        mirrored_df["End"] = pd.to_numeric(end, errors="raise").astype(int)

    merged_df = (
        mirrored_df
        .groupby(["Clone1", "Clone2", "Chromosome"], sort=False, group_keys=False)
        .apply(merge_adjacent_segments)
        .reset_index(drop=True)
    )

    merged_df["length"] = merged_df["End"].astype(int) - merged_df["Start"].astype(int) + 1

    merged_df["size"] = pd.cut(
        merged_df["length"],
        bins=[-1, small_size, mid_size, float("inf")],
        labels=["focal", "medium", "broad"],
    )

    mirrored_subclonal_event_counts = (
        merged_df
        .groupby(["size"])
        .size()
        .reset_index(name="count")
    )

    return mirrored_subclonal_event_counts


def merge_adjacent_segments(df_grp: pd.DataFrame) -> pd.DataFrame:
    df = df_grp.copy()
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)

    # 组内按坐标排序，保证“相邻行”语义正确
    df = df.sort_values(["Start", "End"], kind="mergesort").reset_index(drop=True)

    merged_rows = []
    cur = df.iloc[0].copy()

    for i in range(1, len(df)):
        nxt = df.iloc[i]

        is_contiguous = (cur["End"] == int(nxt["Start"]) - 1)
        same_cna = (cur["Clone1_CNA"] == nxt["Clone1_CNA"])

        if is_contiguous and same_cna:
            # 合并：延长当前段
            cur["End"] = int(nxt["End"])
            # 如果你也希望 region 自动更新，可以最后统一生成；这里先不动
        else:
            merged_rows.append(cur)
            cur = nxt.copy()

    merged_rows.append(cur)

    out = pd.DataFrame(merged_rows)

    return out


def align_tables(bin_df: pd.DataFrame, cn_df: pd.DataFrame):
    """
    Align read depth table (bin_df) and copy number table (cn_df)
    by matching both cell columns and genomic region rows.

    Args:
        bin_df (pd.DataFrame): DataFrame containing read depth per bin.
        cn_df (pd.DataFrame): DataFrame containing predicted copy number per bin.

    Returns:
        tuple: (aligned_bin_df, aligned_cn_df, common_cells, common_regions)
    """

    # ---------- Align columns (cells) ----------
    common_cells = list(set(bin_df.columns[1:]).intersection(cn_df.columns[1:]))
    common_cells.sort()

    if len(common_cells) == 0:
        raise ValueError("No common cell names found. Please check input files.")

    bin_df = bin_df[["region"] + common_cells]
    cn_df = cn_df[["region"] + common_cells]

    # ---------- Align rows (genomic regions) ----------
    common_regions = set(bin_df["region"]).intersection(set(cn_df["region"]))

    bin_df = (
        bin_df[bin_df["region"].isin(common_regions)]
        .sort_values("region")
        .reset_index(drop=True)
    )
    cn_df = (
        cn_df[cn_df["region"].isin(common_regions)]
        .sort_values("region")
        .reset_index(drop=True)
    )

    # Sanity check: ensure identical order
    if not all(bin_df["region"] == cn_df["region"]):
        raise ValueError("Region order mismatch after alignment.")

    print(f"Alignment complete: {len(common_regions)} regions, {len(common_cells)} cells.")

    return bin_df, cn_df, common_cells, common_regions


def compute_rd_cn_l1(bin_path: str, cn_path: str):
    """
    Compute normalized L1 error between read depth and predicted copy number.

    This function:
        1. Reads the two CSV input files.
        2. Aligns them by shared cells and genomic bins.
        3. Parses haplotype copy numbers (hap1|hap2).
        4. Normalizes per cell.
        5. Computes the mean L1 distance between normalized CN and read depth.

    Args:
        bin_path (str): Path to the read depth CSV file (e.g. bin_counts.csv).
        cn_path (str): Path to the predicted copy number CSV file (e.g. haplotype_combined.csv).

    Returns:
        pd.Series: L1 error per bin.
    """

    # ---------- Step 1: Read input files ----------
    bin_df = read_table_auto(bin_path)
    cn_df = read_table_auto(cn_path)

    # ---------- Step 2: Align both tables ----------
    bin_df, cn_df, cells, regions = align_tables(bin_df, cn_df)

    # ---------- Step 3: Parse haplotype copy numbers ----------
    def parse_hap_sum(val):
        """Convert 'hap1|hap2' → hap1 + hap2, return NaN for invalid values."""
        try:
            a, b = val.split("|")
            return float(a) + float(b)
        except Exception:
            return np.nan

    cn_numeric = cn_df[cells].applymap(parse_hap_sum)

    # ---------- Step 4: Convert to numeric matrices ----------
    read_counts = bin_df[cells].astype(float)
    copy_numbers = cn_numeric

    # ---------- Step 5: Normalize each cell by its mean value ----------
    normalized_read_depth = read_counts / read_counts.mean(axis=0)
    normalized_cn = copy_numbers / copy_numbers.mean(axis=0)

    # ---------- Step 6: Compute L1 error per cell ----------
    l1_errors = (normalized_cn - normalized_read_depth).abs().mean(axis=0)
    l1_errors.index = bin_df.columns[1:] 

    print("L1 error calculation completed.")
    return l1_errors

def count_unique_cells(file_path = None):
    df = pd.read_csv(file_path, index_col=0)
    
    unique_cells = df.T.drop_duplicates()
    num_unique_cells = unique_cells.shape[0]

    return num_unique_cells


def long_to_mtx(
    df: pd.DataFrame,
    out_dir: str,
    prefix: str = "cellSNP.tag",
    chr_col: str = "chr",
    pos_col: str = "position",
    cell_col: str = "cell",
    a_col: str = "Acount",
    b_col: str = "Bcount",
    ad_from: str = "B",  # "A" or "B"
    make_vaf: bool = True,
    min_dp = 1,
    min_cells = 1,
):

    os.makedirs(out_dir, exist_ok=True)

    needed = {chr_col, pos_col, cell_col, a_col, b_col}
    miss = needed - set(df.columns)
    if miss:
        raise ValueError(f"df miss: {sorted(miss)}")

    df = df.copy()
    df[pos_col] = df[pos_col].astype(int)
    df[a_col] = df[a_col].astype(float)
    df[b_col] = df[b_col].astype(float)
    df["DP"] = df[a_col] + df[b_col]
    df["AD"] = df[a_col] if ad_from == "A" else df[b_col]

    df = df[df["DP"] > 0].copy()

    df["variant_key"] = df[chr_col].astype(str) + ":" + df[pos_col].astype(str)

    variant_keys = df["variant_key"].unique()
    cell_ids = df[cell_col].unique()

    variant_keys = np.array(sorted(variant_keys))
    cell_ids = np.array(sorted(cell_ids))

    v2i = {v: i for i, v in enumerate(variant_keys)}
    c2j = {c: j for j, c in enumerate(cell_ids)}

    row = df["variant_key"].map(v2i).to_numpy()
    col = df[cell_col].map(c2j).to_numpy()

    AD = coo_matrix((df["AD"].to_numpy(), (row, col)),
                    shape=(len(variant_keys), len(cell_ids))).tocsr()
    DP = coo_matrix((df["DP"].to_numpy(), (row, col)),
                    shape=(len(variant_keys), len(cell_ids))).tocsr()

    ad_path = os.path.join(out_dir, f"{prefix}.AD.mtx")
    dp_path = os.path.join(out_dir, f"{prefix}.DP.mtx")
    scipy.io.mmwrite(ad_path, AD)
    scipy.io.mmwrite(dp_path, DP)

    samples_path = os.path.join(out_dir, f"{prefix}.samples.tsv")
    pd.Series(cell_ids).to_csv(samples_path, sep="\t", header=False, index=False)

    parts = pd.Series(variant_keys).str.split(":", expand=True)
    variants_df = pd.DataFrame({
        "chr": parts[0].values,
        "position": parts[1].astype(int).values,
    })
    variants_path = os.path.join(out_dir, f"{prefix}.variants.tsv")
    variants_df.to_csv(variants_path, sep="\t", index=False)

    vaf_path = f"{out_dir}/{prefix}.VAF.mtx"
    if make_vaf:
        process_variant_data(out_dir, out_dir, min_dp=min_dp, min_cells=min_cells)

    return {
        "AD_path": ad_path,
        "DP_path": dp_path,
        "samples_path": samples_path,
        "variants_path": variants_path,
        "n_variants": int(len(variant_keys)),
        "n_cells": int(len(cell_ids)),
        "AD_nnz": int(AD.nnz),
        "DP_nnz": int(DP.nnz),
    }


def run_cmd(cmd, check=True):
    """Run a shell command and return stdout as text."""
    p = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False
    )
    if check and p.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            f"{' '.join(cmd)}\n\n"
            f"STDERR:\n{p.stderr}"
        )
    return p.stdout


def parse_ad(ad_str):
    """
    Parse FORMAT/AD field exported by bcftools query as a string.
    Expected forms:
      "12,3"   -> (ref=12.0, alt=3.0)
      "." or ".,." or "" -> (None, None)
    """
    if ad_str is None:
        return None, None
    s = str(ad_str).strip()
    if s == "" or s == ".":
        return None, None
    if "," not in s:
        return None, None
    a, b = s.split(",", 1)
    a = a.strip()
    b = b.strip()
    if a == "." or b == "." or a == "" or b == "":
        return None, None
    try:
        return float(a), float(b)
    except ValueError:
        return None, None


def get_samples_from_vcf(vcf_path):
    """Get sample names (cell IDs) from a VCF via bcftools."""
    out = run_cmd(["bcftools", "query", "-l", vcf_path])
    samples = [x.strip() for x in out.splitlines() if x.strip() != ""]
    if len(samples) == 0:
        raise ValueError(f"No samples found in VCF header: {vcf_path}")
    return samples


def check_samples_consistent(vcf_paths, samples0):
    """Ensure all VCFs have the same samples in the same order."""
    for p in vcf_paths[1:]:
        s = get_samples_from_vcf(p)
        if s != samples0:
            raise ValueError(
                "Sample list/order differs across VCFs.\n"
                f"First VCF: {vcf_paths[0]} (n={len(samples0)})\n"
                f"Current VCF: {p} (n={len(s)})\n"
                "You must make sample order consistent before building matrices."
            )


def vcf_chr_files_to_ad_dp_mtx(
    vcf_paths,
    out_dir,
    prefix="cellSNP.tag",
    keep_only_biallelic=True,
    make_vaf: bool = True,
    min_dp = 1,
    min_cells = 1,
):
    """
    Build AD.mtx (alt counts) and DP.mtx (ref+alt) from multiple VCFs (chr1..chr22).
    Uses bcftools query:
      %CHROM %POS %REF %ALT [ %AD ]
    """
    os.makedirs(out_dir, exist_ok=True)

    # 1) samples list from first VCF
    samples = get_samples_from_vcf(vcf_paths[0])
    check_samples_consistent(vcf_paths, samples)
    n_cells = len(samples)

    # 2) COO accumulation
    ad_rows, ad_cols, ad_data = [], [], []
    dp_rows, dp_cols, dp_data = [], [], []
    variants = []  # list of (chr, pos, ref, alt)

    # 3) iterate VCFs
    for vcf_path in vcf_paths:
        fmt = r"%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n"
        out = run_cmd(["bcftools", "query", "-f", fmt, vcf_path])

        for line in out.splitlines():
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue

            chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]

            if keep_only_biallelic and ("," in alt):
                # skip multiallelic records: AD may have >2 entries
                continue

            try:
                pos_i = int(pos)
            except ValueError:
                continue

            ad_fields = parts[4:]
            if len(ad_fields) != n_cells:
                raise ValueError(
                    f"AD columns mismatch in {vcf_path}\n"
                    f"Expected {n_cells} AD fields, got {len(ad_fields)}\n"
                    f"Example line: {line[:200]}..."
                )

            variants.append((chrom, pos_i))
            v_idx = len(variants) - 1

            for j, ad_str in enumerate(ad_fields):
                ref_c, alt_c = parse_ad(ad_str)
                if ref_c is None:
                    continue

                dp = ref_c + alt_c
                if dp <= 0:
                    continue

                # DP always recorded when dp>0
                dp_rows.append(v_idx)
                dp_cols.append(j)
                dp_data.append(dp)

                # AD is alt counts; record only if >0 to keep sparse
                if alt_c > 0:
                    ad_rows.append(v_idx)
                    ad_cols.append(j)
                    ad_data.append(alt_c)

    n_vars = len(variants)

    # 4) build sparse matrices
    AD = coo_matrix((ad_data, (ad_rows, ad_cols)), shape=(n_vars, n_cells)).tocsr()
    DP = coo_matrix((dp_data, (dp_rows, dp_cols)), shape=(n_vars, n_cells)).tocsr()

    # 5) write outputs
    ad_path = os.path.join(out_dir, f"{prefix}.AD.mtx")
    dp_path = os.path.join(out_dir, f"{prefix}.DP.mtx")
    scipy.io.mmwrite(ad_path, AD)
    scipy.io.mmwrite(dp_path, DP)

    samples_path = os.path.join(out_dir, "cellSNP.samples.tsv")
    pd.Series(samples).to_csv(samples_path, sep="\t", header=False, index=False)

    variants_path = os.path.join(out_dir, "cellSNP.variants.tsv")
    pd.DataFrame(variants, columns=["chr", "position"]).to_csv(
        variants_path, sep="\t", index=False
    )

    if make_vaf:
        process_variant_data(out_dir, out_dir, min_dp=min_dp, min_cells=min_cells)


    return {
        "AD_path": ad_path,
        "DP_path": dp_path,
        "samples_path": samples_path,
        "variants_path": variants_path,
        "n_variants": int(n_vars),
        "n_cells": int(n_cells),
        "AD_nnz": int(AD.nnz),
        "DP_nnz": int(DP.nnz),
    }

def intersect_cells_from_cna(tool_hap1_files: List[str], tool_hap2_files: List[str]) -> List[str]:
    """
    Compute the intersection of cell IDs across all tools based on CNA matrices.

    For each tool, we intersect hap1 and hap2 CNA columns (cells), then intersect
    across all tools. The final list is sorted for deterministic column order.
    """
    common = None
    for f1, f2 in zip(tool_hap1_files, tool_hap2_files):
        cnaA = pd.read_csv(f1, index_col=0)
        cnaB = pd.read_csv(f2, index_col=0)
        cells = set(cnaA.columns).intersection(cnaB.columns)
        common = cells if common is None else common.intersection(cells)
    return sorted(list(common)) if common else []


def phase_to_binary(hap1_df: pd.DataFrame, hap2_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    
    h1 = (hap1_df > hap2_df).astype(int)
    h1 = h1.mask(hap1_df == hap2_df, -1).mask((hap1_df.isna() | hap2_df.isna()), -1)
    h2 = (hap2_df > hap1_df).astype(int)
    h2 = h2.mask(hap1_df == hap2_df, -1).mask((hap1_df.isna() | hap2_df.isna()), -1)
    return h1, h2


# def eval_mismatch_switch_gt(g1: pd.DataFrame, g2: pd.DataFrame,
#                             h1: pd.DataFrame, h2: pd.DataFrame, tool_name: str) -> List[Dict]:
#     # g1 = g1.dropna(axis=1, how='all'); g2 = g2.dropna(axis=1, how='all')

#     # idx = g1.index.intersection(h1.index)
#     # col = g1.columns.intersection(h1.columns)
#     # g1, g2, h1, h2 = g1.loc[idx, col], g2.loc[idx, col], h1.loc[idx, col], h2.loc[idx, col]

#     # col = g1.columns

#     # prefixes = col.str.extract(r"(^[^_]+)_")[0].unique()
#     results = []


#     grp_cols = g1.columns.intersection(h1.columns)
#     g1_g, g2_g = g1[grp_cols], g2[grp_cols]
#     h1_g,  h2_g = h1[grp_cols],  h2[grp_cols]

#     # mismatch
#     # mask_1 = (g1_g != -1) & (h1_g != -1)
#     # mask_2 = (g2_g != -1) & (h2_g != -1)

#     mask_1 = g1_g != -1; mask_2 = g2_g != -1

#     print(f"mismatch total diff {(mask_1 != mask_2).sum().sum()}")
#     print(f"mask1: { mask_1.sum().sum()}")
#     print(f"mask2: { mask_2.sum().sum()}")


#     mm = ((h1_g != g1_g) & mask_1).sum().sum()
#     pm = ((h2_g != g2_g) & mask_2).sum().sum()

#     total_m = mask_1.sum().sum()
#     print(f"mask1: { mask_1.sum().sum()}")
#     print(f"mask2: { mask_2.sum().sum()}")


#     # switch error
#     nrow, ncol = g1_g.shape
#     total_sw = 0; sw_1 = 0; sw_2 = 0


#     g1_n = g1_g.to_numpy()
#     g2_n = g2_g.to_numpy()
#     h1_n = h1_g.to_numpy()
#     h2_n = h2_g.to_numpy()

#     g1_cur, g1_next = g1_n[:, :-1], g1_n[:, 1:]
#     g2_cur, g2_next = g2_n[:, :-1], g2_n[:, 1:]
#     h1_cur, h1_next = h1_n[:, :-1], h1_n[:, 1:]
#     h2_cur, h2_next = h2_n[:, :-1], h2_n[:, 1:]

#     valid1 = (g1_cur != -1) & (g1_next != -1)
#     # hap2 相邻 bin 同时有效
#     valid2 = (g2_cur != -1) & (g2_next != -1)


#     sw1_mask = ((g1_cur != h1_cur) | (g1_next != h1_next)) & valid1
#     sw2_mask = ((g2_cur != h2_cur) | (g2_next != h2_next)) & valid2

#     total_sw = valid1.sum()
#     sw_1 = sw1_mask.sum()

#     results.append({
#         "tool_name": tool_name,
#         "mismatch_count": int(mm),
#         "total": int(total_m),
#         "mismatch_ratio": (mm / total_m) if total_m > 0 else None,
#         "switch_error_count": int(sw_1),
#         "total_switch_compare_count": int(total_sw),
#         "switch_error_ratio": (sw_1 / total_sw) if total_sw > 0 else None,
#     })
#     return results




def eval_mismatch_switch_homorozygous_included_predict(g1: pd.DataFrame, g2: pd.DataFrame,
                            h1: pd.DataFrame, h2: pd.DataFrame, tool_name: str) -> List[Dict]:

    results = []


    grp_cols = g1.columns.intersection(h1.columns)
    g1_g, g2_g = g1[grp_cols], g2[grp_cols]
    h1_g,  h2_g = h1[grp_cols],  h2[grp_cols]

    mask_1 = g1_g != -1; mask_2 = g2_g != -1

    print(f"mismatch total diff {(mask_1 != mask_2).sum(axis = 0).shape}")
    print(f"mask1: { mask_1.sum(axis = 0).shape}")
    print(f"mask2: { mask_2.sum(axis = 0).shape}")


    mm = ((h1_g != g1_g) & mask_1).sum(axis = 0)

    total_m = mask_1.sum(axis = 0)


    # switch error
    nrow, ncol = g1_g.shape
    total_sw = 0; sw_1 = 0; sw_2 = 0


    g1_n = g1_g.to_numpy()
    g2_n = g2_g.to_numpy()
    h1_n = h1_g.to_numpy()
    h2_n = h2_g.to_numpy()

    g1_cur, g1_next = g1_n[:-1, :], g1_n[:-1, :]
    g2_cur, g2_next = g2_n[:-1, :], g2_n[:-1, :]
    h1_cur, h1_next = h1_n[:-1, :], h1_n[:-1, :]
    h2_cur, h2_next = h2_n[:-1, :], h2_n[:-1, :]

    valid1 = (g1_cur != -1) & (g1_next != -1)

    gt_change = (g1_cur != g1_next)
    pr_change = (h1_cur != h1_next)

    switches = (gt_change != pr_change) & valid1

    total_sw = valid1.sum(axis=0)
    sw_1 = switches.sum(axis=0)

    results.append({
        "tool_name": tool_name,
        "mismatch_count": int(mm),
        "total": int(total_m),
        "mismatch_ratio": (mm / total_m) if total_m > 0 else None,
        "switch_error_count": int(sw_1),
        "total_switch_compare_count": int(total_sw),
        "switch_error_ratio": (sw_1 / total_sw) if total_sw > 0 else None,
    })
    return results

def eval_mismatch_switch_homorozygous_included(g1: pd.DataFrame, g2: pd.DataFrame,
                            h1: pd.DataFrame, h2: pd.DataFrame, tool_name: str) -> List[Dict]:
    # g1 = g1.dropna(axis=1, how='all'); g2 = g2.dropna(axis=1, how='all')

    # idx = g1.index.intersection(h1.index)
    # col = g1.columns.intersection(h1.columns)
    # g1, g2, h1, h2 = g1.loc[idx, col], g2.loc[idx, col], h1.loc[idx, col], h2.loc[idx, col]

    # col = g1.columns

    # prefixes = col.str.extract(r"(^[^_]+)_")[0].unique()
    results = []


    grp_cols = g1.columns.intersection(h1.columns)
    g1_g, g2_g = g1[grp_cols], g2[grp_cols]
    h1_g,  h2_g = h1[grp_cols],  h2[grp_cols]


    mm = ((h1_g != g1_g)).sum().sum()
    pm = ((h2_g != g2_g)).sum().sum()

    total_m = h1_g.size
    print(f"mask1: { total_m}")
    # print(f"mask2: { mask_2.sum().sum()}")


    # switch error
    nrow, ncol = g1_g.shape
    total_sw = 0; sw_1 = 0; sw_2 = 0


    g1_n = g1_g.to_numpy()
    g2_n = g2_g.to_numpy()
    h1_n = h1_g.to_numpy()
    h2_n = h2_g.to_numpy()

    g1_cur, g1_next = g1_n[:, :-1], g1_n[:, 1:]
    g2_cur, g2_next = g2_n[:, :-1], g2_n[:, 1:]
    h1_cur, h1_next = h1_n[:, :-1], h1_n[:, 1:]
    h2_cur, h2_next = h2_n[:, :-1], h2_n[:, 1:]

    # valid1 = (g1_cur != -1) & (g1_next != -1)
    # # hap2 相邻 bin 同时有效
    # valid2 = (g2_cur != -1) & (g2_next != -1)


    sw1_mask = ((g1_cur != h1_cur) | (g1_next != h1_next)) 
    sw2_mask = ((g2_cur != h2_cur) | (g2_next != h2_next)) 

    total_sw = g2_cur.size
    sw_1 = sw1_mask.sum()

    results.append({
        "tool_name": tool_name,
        "mismatch_count": int(mm),
        "total": int(total_m),
        "mismatch_ratio": (mm / total_m) if total_m > 0 else None,
        "switch_error_count": int(sw_1),
        "total_switch_compare_count": int(total_sw),
        "switch_error_ratio": (sw_1 / total_sw) if total_sw > 0 else None,
    })
    return results


def eval_mismatch_switch_both(g1: pd.DataFrame, g2: pd.DataFrame,
                            h1: pd.DataFrame, h2: pd.DataFrame, tool_name: str) -> List[Dict]:
    # g1 = g1.dropna(axis=1, how='all'); g2 = g2.dropna(axis=1, how='all')

    # idx = g1.index.intersection(h1.index)
    # col = g1.columns.intersection(h1.columns)
    # g1, g2, h1, h2 = g1.loc[idx, col], g2.loc[idx, col], h1.loc[idx, col], h2.loc[idx, col]

    # col = g1.columns

    # prefixes = col.str.extract(r"(^[^_]+)_")[0].unique()
    results = []

    grp_cols = g1.columns.intersection(h1.columns)
    g1_g, g2_g = g1[grp_cols], g2[grp_cols]
    h1_g,  h2_g = h1[grp_cols],  h2[grp_cols]

    # mismatch
    mask_1 = (g1_g != -1) & (h1_g != -1)
    mask_2 = (g2_g != -1) & (h2_g != -1)

    print(f"mismatch total diff {(mask_1 != mask_2).sum().sum()}")
    print(f"mask1: { mask_1.sum().sum()}")
    print(f"mask2: { mask_2.sum().sum()}")


    mm = ((h1_g != g1_g) & mask_1).sum().sum()
    pm = ((h2_g != g2_g) & mask_2).sum().sum()

    total_m = mask_1.sum().sum()
    print(f"mask1: { mask_1.sum().sum()}")
    print(f"mask2: { mask_2.sum().sum()}")


    # switch error
    nrow, ncol = g1_g.shape
    total_sw = 0; sw_1 = 0; sw_2 = 0


    g1_n = g1_g.to_numpy()
    g2_n = g2_g.to_numpy()
    h1_n = h1_g.to_numpy()
    h2_n = h2_g.to_numpy()

    g1_cur, g1_next = g1_n[:, :-1], g1_n[:, 1:]
    g2_cur, g2_next = g2_n[:, :-1], g2_n[:, 1:]
    h1_cur, h1_next = h1_n[:, :-1], h1_n[:, 1:]
    h2_cur, h2_next = h2_n[:, :-1], h2_n[:, 1:]

    valid1 = (g1_cur != -1) & (g1_next != -1) & (h1_cur != -1) & (h1_next != -1)
    # hap2 相邻 bin 同时有效
    valid2 = (g2_cur != -1) & (g2_next != -1) & (h2_cur != -1) & (h2_next != -1)


    sw1_mask = ((g1_cur != h1_cur) | (g1_next != h1_next)) & valid1
    sw2_mask = ((g2_cur != h2_cur) | (g2_next != h2_next)) & valid2

    total_sw = valid1.sum()
    sw_1 = sw1_mask.sum()

    results.append({
        "tool_name": tool_name,
        "mismatch_count": int(mm),
        "total": int(total_m),
        "mismatch_ratio": (mm / total_m) if total_m > 0 else None,
        "switch_error_count": int(sw_1),
        "total_switch_compare_count": int(total_sw),
        "switch_error_ratio": (sw_1 / total_sw) if total_sw > 0 else None,
    })
    return results