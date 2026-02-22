from collections import Counter
import json
import os
import pandas as pd
import numpy as np
import scipy
from sklearn.metrics import adjusted_rand_score as adjustedRandIndex, mean_squared_error
from sklearn.metrics import adjusted_mutual_info_score as AMI
from scipy.sparse import csr_matrix, diags


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
    bin_df = pd.read_csv(bin_path)
    cn_df = pd.read_csv(cn_path)

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

def get_ref_path(filename):
    return os.path.join(os.path.dirname(os.path.dirname(__file__)), "ref", filename)

def load_centromeres(genome="hg38"):
    print(f"Loading {genome} centromere data")
    df = pd.read_csv(get_ref_path(f"{genome}_centromeres.csv"))
    return df

def annotate_segments(df, small_size: int = 5_000_000,mid_size: int = 20_000_000,
        genome: str = "hg38",):
    
    centromeres_df = load_centromeres(genome)
    records = []

    for region, row in df.iterrows():
        chrom, coords = region.split(":")
        start, end = map(int, coords.split("-"))

        for cell, val in row.items():
            if pd.isna(val):
                continue
            sval = str(val)
            if sval in ("1|1", "nan|nan"):
                continue

            size_label = classify_size(start, end, small_size, mid_size)
            type_label = classify_whole_arm_events(
                sval, chrom, start, end, centromeres_df
            )
            records.append(
                {
                    "cell": cell,
                    "region": region,
                    "value": sval,
                    "size": size_label,
                    "type": type_label,
                }
            )

    return pd.DataFrame(records)

def classify_size(start, end, small_size, mid_size):
        size = int(end) - int(start) + 1
        if size <= small_size:
            return "small"
        elif size <= mid_size:
            return "middle"
        return "large"


def classify_cnv(val: str) -> str:
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

def classify_whole_arm_events(
    val: str, chrom: str, start: int, end: int, centromeres_df: pd.DataFrame
) -> str:
    try:
        left, right = map(int, str(val).split("|"))
    except Exception:
        return "UNKNOWN"

    ch = str(chrom)
    if not ch.startswith("chr"):
        ch = "chr" + ch

    row = centromeres_df.loc[centromeres_df["chrom"] == ch]
    if row.empty:
        return classify_cnv(val)
    row = row.iloc[0]

    p_end = int(row["p_arm_end"])
    q_start = int(row["q_arm_start"])
    chr_len = int(row["chrom_length"])

    covers_p = (start <= 1) and (end >= p_end)
    covers_q = (start <= q_start) and (end >= chr_len)

    if (left == 0 or right == 0) and (covers_p or covers_q):
        return "WCL"  # Whole-arm loss
    if (left > 1 or right > 1) and (covers_p or covers_q):
        return "WGD"  # Whole-arm gain

    return classify_cnv(val)


def get_top_clusters(input_file, n=5, small=True):

    cluster_col = "clone_id"

    data = pd.read_csv(input_file)

    clusters = data[cluster_col].astype(str)  
    cluster_counts = clusters.value_counts()

    if small:
        result = cluster_counts.nsmallest(n)
    else:
        result = cluster_counts.nlargest(n)


    return result


def evaluate_clustering_results(path1, path2, tool1, tool2):

    data1 = pd.read_csv(path1)

    data2 = pd.read_csv(path2)

    # Merge data on 'cell_id'
    data = pd.merge(data1, data2, on='cell_id', suffixes=('_1', '_2'))

    data = data.dropna(subset=["clone_id_1", "clone_id_2"])

    
    clusters1 = data['clone_id_1']

    clusters2 = data['clone_id_2']

    
    # Compute the Adjusted Rand Index (ARI) between clusters and cell clones
    ari = adjustedRandIndex(clusters1, clusters2)
    
    # Compute the Adjusted Mutual Information (AMI) between clusters and cell clones
    ami = AMI(clusters1, clusters2)
    
    # Return the results as a dictionary
    return ari, ami

def count_unique_cells(file_path = None):
    df = pd.read_csv(file_path, index_col=0)
    
    unique_cells = df.T.drop_duplicates()
    num_unique_cells = unique_cells.shape[0]

    return num_unique_cells

def split_region(region_str):
    try:
        chrom, coords = region_str.split(":")
        start, end = map(int, coords.split("-"))
        return chrom, start, end
    except Exception as e:
        raise ValueError(f"Invalid segment format: {region_str}") from e


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
        [master[["chrom","start","end","region"]],
        pd.DataFrame(np.nan, index=master.index, columns=value_cols)],
        axis=1
    )

    for chrom in master.chrom.unique():
        m = master[master.chrom == chrom]
        o = old_df[old_df.chrom == chrom]

        for idx, row in m.iterrows():
            s1, e1 = row.start, row.end

            hits = o[(o.start <= s1) & (o.end >= e1)]
            if len(hits) !=1:
                continue

            hit = hits.iloc[0]

            out.loc[idx, value_cols] = hit[value_cols].values

    return out


def align_cna_bins(cna_df_1,cna_df_2):
    cna_df1 = cna_df_1.copy()
    cna_df2 = cna_df_2.copy()


    cna_df1[['chrom', 'start', 'end']] = cna_df1['region'].apply(
        lambda s: pd.Series(split_region(s))
    )


    cna_df2[['chrom', 'start', 'end']] = cna_df2['region'].apply(
        lambda s: pd.Series(split_region(s))
    )

    df1 = cna_df1[['chrom', 'start', 'end']]
    df2 = cna_df2[['chrom', 'start', 'end']] 

    master = build_master_bins(df1, df2)
    master["region"] = master.apply(
        lambda r: f"{r.chrom}:{r.start}-{r.end}", axis=1
    )

    cna_df1_new = map_to_master(master,cna_df1)
    cna_df1_new.drop(columns=['chrom', 'start', 'end'],inplace=True)

    cna_df2_new = map_to_master(master,cna_df2)
    cna_df2_new.drop(columns=['chrom', 'start', 'end'],inplace=True)


    return cna_df1_new, cna_df2_new


def create_lazac_input(path,out_path):
    # node,chrom,start,end,cn_a,cn_b
    cna_df = pd.read_csv(path,index_col = 0)

    cna_df = cna_df.T

    cna_df.index.name = "node"

    long_df = cna_df.reset_index().melt(
        id_vars = 'node',
        var_name = "region",
        value_name = "cn"
    )
    
    long_df[['chrom', 'start', 'end']] = long_df['region'].apply(
        lambda s: pd.Series(split_region(s))
    )

    long_df['cn_a'] = long_df['cn'].str.split("|").str[0]
    long_df['cn_b'] = long_df['cn'].str.split("|").str[1]

    long_df.drop(columns=["region",'cn'],inplace=True)

    empty_like = {"", "NA", "<NA>", "None", "nan"}
    long_df["cn_a"] = long_df["cn_a"].replace(list(empty_like), pd.NA)
    long_df["cn_b"] = long_df["cn_b"].replace(list(empty_like), pd.NA)

    long_df = long_df.dropna(subset=["cn_a", "cn_b"])


    long_df.to_csv(out_path,index=False)

    return 


def get_final_parsimony_score(json_file):
    with open(json_file, "r") as f:
        data = json.load(f)

    # last_iter = data[-1]
    # scores = last_iter["scores"]

    
    # counter = Counter(scores)
    final_score =  min(s
        for it in data
        for s in it["scores"]
        if s >= 0
    )

    return final_score


def map_region_to_variants(variant_info,cna_df):
    """
    Returns a dictionary mapping genomic regions to the indices of variants located within them.
    Example:
    {
        "chr1:1000-2000": [0, 1, 5],   # Indices of rows in variant_info inside this region
        "chr2:5000-6000": [10, 12]
    }

    """

    region_to_var_rows = {}
    print("Mapping regions to variants...")

    unique_chroms = variant_info['chrom'].unique()

    for chrom in unique_chroms:
        chrom_vars = variant_info[variant_info['chrom'] == str(chrom)]

        if chrom_vars.empty: continue
        
        var_positions = chrom_vars['pos'].values
        var_orig_indices = chrom_vars.index.values 

        for region in cna_df.index:

            r_chrom, r_start, r_end = split_region(region)

            if r_chrom == str(chrom):

                s = np.searchsorted(var_positions, r_start)
                e = np.searchsorted(var_positions, r_end)
                
                if s < e:
                    region_to_var_rows[region] = var_orig_indices[s:e]


    print(f"Mapped {len(region_to_var_rows)} regions containing variants.")

    return region_to_var_rows




def extract_vaf_by_binary_mask(target_hap, cna_df, VAF, region_to_var_rows, cell_list):

    mask_df = (cna_df == target_hap)
        
    cell_name_to_col_idx = {name: i for i, name in enumerate(cell_list)}

    collected_values = []

    for region_name, var_rows in region_to_var_rows.items():
        if region_name not in mask_df.index:
                continue
        
        region_mask = mask_df.loc[region_name]

        target_cell_names = region_mask[region_mask].index


        if len(target_cell_names) == 0:
                continue
        
        target_col_indices = [cell_name_to_col_idx[c] for c in target_cell_names if c in cell_name_to_col_idx]

        if not target_col_indices:
                continue
        
        sub_matrix = VAF[var_rows, :][:, target_col_indices]

        data = sub_matrix.data

        if len(data) > 0:
            collected_values.append(data)

    if len(collected_values) == 0:
        return np.array([])

    return np.concatenate(collected_values)

def sparse_divide(A, B):

    B_inv = B.copy()
    B_inv.data = 1.0 / B_inv.data
    return A.multiply(B_inv)

def process_variant_data(base_input_path, output_dir, min_dp=2, min_cells=3, prefix = "cellSNP"):

    os.makedirs(output_dir, exist_ok=True)
    
    print("-" * 30)
    print(f"Processing data from: {base_input_path}")
    print(f"Output directory: {output_dir}")


    print("Loading AD and DP matrices...")
    ad_path = os.path.join(base_input_path, f"{prefix}.AD.mtx")
    dp_path = os.path.join(base_input_path, f"{prefix}.DP.mtx")
    
    AD = scipy.io.mmread(ad_path).tocsr()
    DP = scipy.io.mmread(dp_path).tocsr()
    print(f"Matrix shape: {AD.shape}")


    VAF = sparse_divide(AD, DP)
    print(f"Initial Non-zero elements in VAF: {VAF.nnz}")

    print(f"Applying filter: Masking elements where DP < {min_dp} ...")
    dp_mask = (DP >= min_dp)

    VAF = VAF.multiply(dp_mask)
    print(f"Non-zero elements in VAF after DP filter: {VAF.nnz}")


    print(f"Applying filter: Zeroing out variants appearing in < {min_cells} cells...")
    dp_nonzero = (DP > 0)
    cells_with_coverage_count = np.array(dp_nonzero.sum(axis=1)).flatten()

    rows_to_keep_mask = (cells_with_coverage_count >= min_cells).astype(int)
    diag_filter = diags(rows_to_keep_mask)

    AD = AD.multiply(dp_mask)
    DP = DP.multiply(dp_mask)



    VAF = diag_filter.dot(VAF).tocsr()
    AD = diag_filter.dot(AD).tocsr()
    DP = diag_filter.dot(DP).tocsr()
    print(f"Non-zero elements in VAF after cell filter: {VAF.nnz}")



    print("Writing filtered matrices and samples to disk...")


    scipy.io.mmwrite(os.path.join(output_dir, f"{prefix}.VAF.filtered.mtx"), VAF)
    scipy.io.mmwrite(os.path.join(output_dir, f"{prefix}.AD.filtered.mtx"), AD)
    scipy.io.mmwrite(os.path.join(output_dir, f"{prefix}.DP.filtered.mtx"), DP)
    

    print("Processing complete.")
    print("-" * 30)


def match_snvs_to_bins(snv_df, cna_df):
    """
    将 SNV 映射到对应的 Bin 索引。
    
    参数:
    snv_df : DataFrame, 包含列 ['chrom', 'pos']
    cna_df : DataFrame, 包含列 ['region']
    
    返回:
    snv_to_bin_idx : np.array, 长度等于 snv_df 的行数。
                     值为对应的 bin 在 cna_df 中的 index。
                     如果没找到对应 bin，值为 -1。
    """

    snv_to_bin_idx = np.full(len(snv_df), -1, dtype=int)
    
    chromosomes = snv_df['chrom'].unique()

    cna_df[['chrom', 'start', 'end']] = cna_df['region'].apply(
        lambda s: pd.Series(split_region(s))
    )
    
    for chrom in chromosomes:

        s_mask = (snv_df['chrom'] == chrom)
        b_mask = (cna_df['chrom'] == chrom)
        
        if not b_mask.any(): continue

        curr_snvs = snv_df[s_mask]
        curr_bins = cna_df[b_mask]
        
        intervals = pd.IntervalIndex.from_arrays(curr_bins['start'], curr_bins['end'], closed='both')

        relative_indices = intervals.get_indexer(curr_snvs['pos'])

        valid_mask = relative_indices != -1
        
        global_indices = curr_bins.index[relative_indices[valid_mask]].values
        
        snv_global_indices = np.where(s_mask)[0]
        snv_to_bin_idx[snv_global_indices[valid_mask]] = global_indices
        
    return snv_to_bin_idx


def calc_method_ll_from_bins(snv_v, snv_r, snv_indices, cn_A, cn_B,
                             clip_min=0.05, clip_max=0.95):

    N_SNV = snv_v.shape[0]
    ll_per_snv = np.full(N_SNV, np.nan, dtype=float) 

    valid_mask = (snv_indices != -1)
    if not np.any(valid_mask):
        return ll_per_snv  

    v = snv_v[valid_mask]          # (N_valid, N_Cell)
    r = snv_r[valid_mask]          # (N_valid, N_Cell)
    idx = snv_indices[valid_mask]  # (N_valid,)

    
    X_A = cn_A[idx]   # (N_valid, N_Cell)
    X_B = cn_B[idx]

    total_cn = X_A + X_B

    with np.errstate(divide='ignore', invalid='ignore'):

        prob_v_A = X_A / total_cn
        prob_r_A = X_B / total_cn  

        prob_v_A = np.clip(prob_v_A, clip_min, clip_max)
        prob_r_A = np.clip(prob_r_A, clip_min, clip_max)

        ll_A = np.sum(v * np.log(prob_v_A) + r * np.log(prob_r_A), axis=1)


        prob_v_B = prob_r_A
        prob_r_B = prob_v_A

        ll_B = np.sum(v * np.log(prob_v_B) + r * np.log(prob_r_B), axis=1)


    ll_valid = np.maximum(ll_A, ll_B)  # (N_valid,)

    
    ll_per_snv[valid_mask] = ll_valid

    return ll_per_snv

def bootstrap_llr(snv_scores, n_boot=1000, random_state=None):

    snv_scores = np.asarray(snv_scores)
    n = snv_scores.shape[0]

    rng = np.random.default_rng(random_state)

    boot_totals = np.empty(n_boot, dtype=float)

    for b in range(n_boot):
        
        idx = rng.integers(0, n, size=n)
        boot_totals[b] = snv_scores[idx].sum()

    return boot_totals

def permuted_llr_once(
    snv_v, snv_r,
    snv_indices1, snv_indices2,
    cn1_A, cn1_B, cn2_A, cn2_B,
    rng
):

    perm1 = rng.permutation(cn1_A.shape[0])

    cn1_Ap = cn1_A[perm1, :]
    cn1_Bp = cn1_B[perm1, :]

    perm2 = rng.permutation(cn2_A.shape[0])
    cn2_Ap = cn2_A[perm2, :]
    cn2_Bp = cn2_B[perm2, :]

    ll1 = calc_method_ll_from_bins(snv_v, snv_r, snv_indices1, cn1_Ap, cn1_Bp)
    ll2 = calc_method_ll_from_bins(snv_v, snv_r, snv_indices2, cn2_Ap, cn2_Bp)

    mask = (~np.isnan(ll1)) & (~np.isnan(ll2))
    return np.sum(ll1[mask] - ll2[mask])


def permutation_llr(
    snv_v, snv_r,
    snv_indices1, snv_indices2,
    cn1_A, cn1_B, cn2_A, cn2_B,
    n_perm=2000, random_state=0
):
    rng = np.random.default_rng(random_state)
    llr_null = np.zeros(n_perm)

    for i in range(n_perm):
        llr_null[i] = permuted_llr_once(
            snv_v, snv_r,
            snv_indices1, snv_indices2,
            cn1_A, cn1_B, cn2_A, cn2_B,
            rng
        )

    return llr_null

def calculate_llr(
    AD_df1, DP_df1,
    AD_df2, DP_df2,
    snv_indices1, snv_indices2,
    cna1_A_df, cna1_B_df,   
    cna2_A_df, cna2_B_df,    
    n_perm = 2000
):

    idx1 = AD_df1.index
    idx2 = AD_df2.index

    common_snv_index = idx1.intersection(idx2)

    pos1 = idx1.get_indexer(common_snv_index)
    pos2 = idx2.get_indexer(common_snv_index)


    common_cells = (
        AD_df1.columns
        .intersection(DP_df1.columns)
        .intersection(AD_df2.columns)
        .intersection(DP_df2.columns)
        .intersection(cna1_A_df.columns)
        .intersection(cna1_B_df.columns)
        .intersection(cna2_A_df.columns)
        .intersection(cna2_B_df.columns)
    )
    common_cells = list(common_cells)

    if len(common_cells) == 0:
        return {
            "llr_obs": np.nan,
            "p_value": np.nan,
            "boot_totals": [],
            "n_common_cells": 0,
            "n_common_snvs": 0,
            "n_snvs_used": 0,
            "note": "No common cells across inputs",
        }

    snv_v1 = AD_df1.loc[common_snv_index, common_cells].to_numpy()
    snv_r1 = DP_df1.loc[common_snv_index, common_cells].to_numpy() - snv_v1

    cn1_A = cna1_A_df.loc[:,common_cells].to_numpy()
    cn1_B = cna1_B_df.loc[:,common_cells].to_numpy()


    snv_v2 = AD_df2.loc[common_snv_index, common_cells].to_numpy()
    snv_r2 = DP_df2.loc[common_snv_index, common_cells].to_numpy() - snv_v2
    cn2_A = cna2_A_df.loc[:,common_cells].to_numpy()
    cn2_B = cna2_B_df.loc[:,common_cells].to_numpy()

    print("cn2_B shape:", cn2_B.shape)
    print("before snv_indices2 min/max:", int(snv_indices2.min()), int(snv_indices2.max()))
    print("unique bins in snv_indices2:", len(np.unique(snv_indices2)))

    snv_indices1 = snv_indices1[pos1]
    snv_indices2 = snv_indices2[pos2]

    print("cn2_B shape:", cn2_B.shape)
    print("after snv_indices2 min/max:", int(snv_indices2.min()), int(snv_indices2.max()))
    print("unique bins in snv_indices2:", len(np.unique(snv_indices2)))

    ll1 = calc_method_ll_from_bins(snv_v1, snv_r1, snv_indices1, cn1_A, cn1_B)
    ll2 = calc_method_ll_from_bins(snv_v2, snv_r2, snv_indices2, cn2_A, cn2_B)

    # # 需要修改成两个工具共同的snv
    mask = (~np.isnan(ll1)) & (~np.isnan(ll2))
    ll1 = ll1[mask]
    ll2 = ll2[mask]
    snv_scores = ll1 - ll2
    llr_obs = snv_scores.sum()
    print(f"ll1 : {ll1[:5]}")
    print(f"ll2 : {ll2[:5]}")


    boot_totals = bootstrap_llr(snv_scores, n_boot=n_perm, random_state=0)

    if llr_obs >= 0:
       
        extreme = np.sum(boot_totals <= 0)
    else:

        extreme = np.sum(boot_totals >= 0)

    p_val = (extreme + 1) / (len(boot_totals) + 1)  

    return {
        "llr_obs": llr_obs,
        # "boot_mean": boot_totals.mean(),
        # "boot_std": boot_totals.std(ddof=1),
        "p_value": p_val,
        "boot_totals": boot_totals.tolist()
    }



