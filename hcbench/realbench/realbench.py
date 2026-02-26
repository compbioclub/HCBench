# hcbench/realbench/realbench.py
import gzip
from operator import gt
import os
import ast
import numpy as np
import pandas as pd
from typing import List, Tuple, Dict, Optional

import scipy
from .utils import align_cna_bins, align_tables, calculate_llr, compute_rd_cn_l1, count_unique_cells, create_lazac_input, evaluate_clustering_results, extract_vaf_by_binary_mask, get_final_parsimony_score, get_top_clusters, map_region_to_variants, match_snvs_to_bins
from hcbench.gtbench import gtbench
from ..utils import align, annotate_segments, categorize_and_save, check_binsize, eval_mismatch_switch_both, eval_mismatch_switch_gt, eval_mismatch_switch_homorozygous_included, evaluate_haplotype_predictions, get_cell_profile_size, get_cluster_size, get_seg_cna_event_num, get_segment_metric, get_segment_metric_both, get_segment_overlap_ratio, intersect_cells_from_cna, phase_to_binary, process_folder_for_metrics, read_and_drop_empty, real_cell_mismatch_error, real_cell_switch_error
import subprocess
from hcbench.parsers.utils import map_cell_to_barcode
import itertools
from scipy.stats import ttest_rel

class RealBench:

    def __init__(self, output_dir: str = "./output"):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def rddetect(
        self,
        tool_cna_files: List[str],
        bin_count_files: List[str],
        tool_names: List[str],
        outfile: str = "rd_cn_l1.csv",
        ttest_outfile: str = "rd_cn_l1.paired_ttest.csv",
        ttest_alternative: str = "two-sided",   # "two-sided" / "less" / "greater"
    ) -> pd.DataFrame:

        results = {}
        check_binsize(tool_cna_files, tool_names, strict=False)

        for path, bin_count_file,name in zip(tool_cna_files, bin_count_files, tool_names):
            l1_error = compute_rd_cn_l1(bin_count_file, path)
    
            results[name] = l1_error

        df = pd.concat(results, axis=1)
        df.columns = tool_names  # optional, to enforce column order
        df.index.name = "Cells"

        # --- Save output ---
        out = os.path.join(self.output_dir, outfile)
        df.to_csv(out, index=True)  # index=region will be written as first column
        print(f"Region-wise L1 error table saved to {out}")

        common_cells = df.dropna().index
        df = df.loc[common_cells]
        ttest_rows = []
        for a in tool_names:
            for b in tool_names:
                if a == b:
                    continue
                if a not in df.columns or b not in df.columns:
                    continue
                sub = df[[a, b]].dropna()
                if sub.shape[0] < 2:
                    continue

                x = sub[a].to_numpy(dtype=float)
                y = sub[b].to_numpy(dtype=float)

                res = ttest_rel(x, y, alternative=ttest_alternative)

                t_stat = np.asarray(res.statistic).reshape(-1)[0]
                p_val  = np.asarray(res.pvalue).reshape(-1)[0]

                ttest_rows.append({
                    "tool_a": a,
                    "tool_b": b,
                    "n_cells": int(sub.shape[0]),
                    "mean_a": float(np.mean(x)),
                    "mean_b": float(np.mean(y)),
                    "mean_diff_a_minus_b": float(np.mean(x - y)),
                    "t": float(t_stat),
                    "p": float(p_val),
                })

        ttest_df = pd.DataFrame(ttest_rows)
        ttest_out = os.path.join(self.output_dir, ttest_outfile)
        ttest_df.to_csv(ttest_out, index=False)
        print(f"Paired t-test table saved to {ttest_out}")


    def rddetect_segment():
        pass
    
    def cloneSize(
        self,
        tool_clone_files: List[str],
        tool_names: List[str],
        top_n: int = 5,
        is_small = True,
        outfile: str = "clone_size.csv",
    ) -> pd.DataFrame:

        result_list = []

        for path, name in zip(tool_clone_files, tool_names):
            cluster_counts = get_top_clusters(path, top_n, is_small)

            df = pd.DataFrame({
                "Tool": [name] * top_n,
                "Num": cluster_counts.values.astype(int),
                "Rank": range(1, top_n + 1)
                })

            result_list.append(df)

        final_df = pd.concat(result_list, ignore_index=True)

         # --- Save output ---
        out = os.path.join(self.output_dir, outfile)
        final_df.to_csv(out, index=False) 
        print(f"clone size saved to {out}")

        return final_df
    
    def clusterConsistency(
        self,
        tool_cluster_files: List[str],
        tool_names: List[str],
    ):
        result_df_ari = pd.DataFrame(columns=tool_names, index=tool_names)
        result_df_ami = pd.DataFrame(columns=tool_names, index=tool_names)


        for path1, tool1 in zip(tool_cluster_files, tool_names):
            for path2, tool2 in zip(tool_cluster_files, tool_names):

                ari, ami = evaluate_clustering_results(path1, path2, tool1, tool2)
                
                result_df_ari.loc[tool1, tool2] = ari
                result_df_ari.loc[tool2, tool1] = ari

                result_df_ami.loc[tool1, tool2] = ami
                result_df_ami.loc[tool2, tool1] = ami

        out_ari = os.path.join(self.output_dir, "clustering_ARI.csv")
        result_df_ari.to_csv(out_ari)
        print(f"Clustering ARI matrix saved to {out_ari}") 
        out_ami = os.path.join(self.output_dir, "clustering_AMI.csv")
        result_df_ami.to_csv(out_ami)
        print(f"Clustering AMI matrix saved to {out_ami}")

        return result_df_ari, result_df_ami

    def cellprofile(
        self,
        tool_cna_files: List[str],
        tool_names: List[str],
        outfile = "unique_cell_profile.csv"
    ):
        
        results = {}

        for path, name in zip(tool_cna_files, tool_names):
            num_unique_cells = count_unique_cells(path)

            results[name] = num_unique_cells

        df = pd.DataFrame.from_dict(results, orient="index", columns=["unique_cells_profile"])
        df.index.name = "Tool"
        out = os.path.join(self.output_dir, outfile)
        df.to_csv(out)  # index=region will be written as first column
        print(f"Unique Cell Profile is saved to {out}")


    def cndetect(
        self,
        tool_cna_files: List[str],
        tool_names: List[str],
        haplotype: Optional[str] = "combined",
        outfile_prefix: str = "bin_level",
        index_col: Optional[int] = 0,
    ) -> pd.DataFrame:

        result_df_rmse = pd.DataFrame(columns=tool_names, index=tool_names)
        result_df_scc = pd.DataFrame(columns=tool_names, index=tool_names)
        result_df_acc = pd.DataFrame(columns=tool_names, index=tool_names)

        for path1, tool1 in zip(tool_cna_files, tool_names):
            for path2, tool2 in zip(tool_cna_files, tool_names):
                print(f"{tool1}----{tool2}")

                truth = pd.read_csv(path1)
                pred = pd.read_csv(path2)
                print("align before", pred.shape)

                truth,pred = align_cna_bins(truth,pred)

                print("align after", pred.shape)


                rmse, scc, acc = evaluate_haplotype_predictions(pred, truth, haplotype)

                result_df_rmse.loc[tool1, tool2] = rmse
                result_df_scc.loc[tool1, tool2] = scc
                result_df_acc.loc[tool1, tool2] = acc


        out_rmse = os.path.join(self.output_dir, f"{outfile_prefix}_rmse.csv")
        result_df_rmse.to_csv(out_rmse)
        print(f"rmse matrix saved to {out_rmse}") 
        out_scc = os.path.join(self.output_dir, f"{outfile_prefix}_scc.csv")
        result_df_scc.to_csv(out_scc)
        print(f"scc matrix saved to {out_scc}") 
        out_acc = os.path.join(self.output_dir, f"{outfile_prefix}_acc.csv")
        result_df_acc.to_csv(out_acc)
        print(f"acc matrix saved to {out_acc}") 

    def dolazactree( self,
        tool_cna_files: List[str],
        tool_names: List[str],
        outfile : Optional[str] = "parsimony_score.csv"
    ):
        
        results = []

        for path, name in zip(tool_cna_files, tool_names):
            out_tool_dir = os.path.join(self.output_dir, f"lazac_{name}")
            os.makedirs(out_tool_dir, exist_ok=True)

            cn_profile_file = f"{out_tool_dir}/{name}_cn_profile.csv"

            create_lazac_input(path,cn_profile_file)

            cmd = [
                "lazac", 
                "nni", 
                cn_profile_file,
                "-a", "2",
                "-o", f"{out_tool_dir}/{name}"
            ]

            print("Running:", " ".join(cmd))
            subprocess.run(cmd, check=True)
            print("Finished!")

            score = get_final_parsimony_score(f"{out_tool_dir}/{name}_info.json")

            results.append({
                "Tool": name,
                "Parsimony": score
            })

        df = pd.DataFrame(results)
        out = os.path.join(self.output_dir, outfile)
        df.to_csv(out,index=False)  
        print(f"Parsimony Score is saved to {out}")

        return df
    

    def findVAFDistribution(
        self,
        tool_cna_files: List[str],
        tool_names: List[str],
        VAF_path_list : str,
        cell_lists : List[str],
        variant_pos_df_list : pd.DataFrame,
        target_cn_status : str
    ):
        
        
        parts = target_cn_status.split("|")
        target_cn_1 = parts[0]
        target_cn_2 = parts[1] if len(parts) > 1 else parts[0]

        is_symmetric = (target_cn_1 == target_cn_2)

        for path, name, VAF_path,cell_list ,variant_pos_df in zip(tool_cna_files, tool_names, VAF_path_list , cell_lists, variant_pos_df_list):

            print(f"Loading VAF matrix from {VAF_path}...")
            VAF = scipy.io.mmread(VAF_path).tocsr()

            cna_df = pd.read_csv(path,index_col=0)

            region_to_var_rows = map_region_to_variants(variant_pos_df, cna_df)

            vafs_A = extract_vaf_by_binary_mask(target_cn_1 + "|"+ target_cn_2, cna_df, VAF, region_to_var_rows,cell_list)

            total_vafs = vafs_A

            if not is_symmetric:
                vafs_B = extract_vaf_by_binary_mask(f"{target_cn_2}|{target_cn_1}", cna_df, VAF, region_to_var_rows, cell_list)
                total_vafs = np.concatenate([vafs_A, vafs_B])

            output_filename = f"{name}_vafs_{target_cn_1}_{target_cn_2}.csv"
            save_path = os.path.join(self.output_dir, output_filename)

            pd.DataFrame(total_vafs, columns=["VAF"]).to_csv(save_path, index=False)
            print(f"  -> Saved to {save_path}")


    @staticmethod
    def load_vcf_pos(vcf_path):
        print(f"Loading VCF from {vcf_path}...")
        positions = []
        chroms = []
        
        if vcf_path.endswith(".vcf.gz"):
            with gzip.open(vcf_path, 'rt') as f:
                for line in f:
                    if line.startswith('#'): continue
                    parts = line.strip().split('\t')
                    chroms.append(parts[0]) # CHROM
                    positions.append(int(parts[1])) # POS
        elif vcf_path.endswith(".tsv"):
            pass
                
        return pd.DataFrame({'chrom': chroms, 'pos': positions})
    
    

    def calLLR(
        self,
        tool_hap1_files: List[str],
        tool_hap2_files: List[str],
        tool_names: List[str],
        snv_paths: List[str],     
        cell_lists: List[List[str]],
        variant_pos_dfs: List[pd.DataFrame]
    ):
        # 1) common cells 只算一次
        common_cells = intersect_cells_from_cna(tool_hap1_files, tool_hap2_files)
        if len(common_cells) == 0:
            raise ValueError("The intersection of cells across all tools is empty.")
        print(f"Common cells across all tools: {len(common_cells)}")


        results = []
        cnaA_list = []
        cnaB_list = []
        snv_to_bin_list = []
        AD_list = []
        DP_list = []

        # 2) Outer loop: treat each tool as the GT (determines which SNV matrix to use)
        for gt_idx, gt_name in enumerate(tool_names):
            gt_snv_path = snv_paths[gt_idx]
            gt_cell_list = cell_lists[gt_idx]
            gt_variant_df = variant_pos_dfs[gt_idx]

            cnaA = pd.read_csv(tool_hap1_files[gt_idx], index_col=0)
            cnaB = pd.read_csv(tool_hap2_files[gt_idx], index_col=0)

            print(f"\n=== GT tool: {gt_name} ===")

            # 2.1 Load GT's AD/DP sparse matrices
            AD = scipy.io.mmread(f"{gt_snv_path}/cellSNP.AD.filtered.mtx").tocsr()
            DP = scipy.io.mmread(f"{gt_snv_path}/cellSNP.DP.filtered.mtx").tocsr()

            # Filter out variants that are all-zero across cells (both AD and DP)
            non_empty_mask = (AD.getnnz(axis=1) > 0) | (DP.getnnz(axis=1) > 0)

            # Sanity check: variant_pos_df must align with unfiltered AD/DP rows
            if len(gt_variant_df) != AD.shape[0]:
                raise ValueError(
                    f"[{gt_name}] variant_pos_df length ({len(gt_variant_df)}) does not match "
                    f"AD/DP row count ({AD.shape[0]}). Ensure they are aligned."
                )

            AD_f = AD[non_empty_mask]
            DP_f = DP[non_empty_mask]
            variant_f = gt_variant_df.loc[non_empty_mask].reset_index(drop=True)

            print(f"[{gt_name}] Total variants: {AD.shape[0]}, kept: {AD_f.shape[0]}")

            variant_key = (
                variant_f["chrom"].astype(str) + ":" +
                variant_f["pos"].astype(str)
            )


            # 2.2 Build AD_df/DP_df with GT's column order, then subset to common cells
            AD_df_all = pd.DataFrame.sparse.from_spmatrix(
                AD_f, index=variant_key, columns=gt_cell_list
            )
            DP_df_all = pd.DataFrame.sparse.from_spmatrix(
                DP_f, index=variant_key, columns=gt_cell_list
            )

            missing_in_gt = [c for c in common_cells if c not in AD_df_all.columns]
            if missing_in_gt:
                raise ValueError(
                    f"[{gt_name}] GT SNV matrices are missing {len(missing_in_gt)} common cells "
                    f"(e.g., {missing_in_gt[:10]}). Ensure SNV and CNA cell IDs are consistent."
                )

            # Restrict to the shared cells and enforce deterministic order
            AD_df = AD_df_all.loc[:, common_cells]
            DP_df = DP_df_all.loc[:, common_cells]

            if not AD_df.index.is_unique:
                dup_n = int(AD_df.index.duplicated(keep=False).sum())
                print(f"[{gt_name}] Found duplicated variant_key rows: {dup_n}. Summing duplicates...")

                # (a) AD/DP 对重复 key 求和
                AD_df = AD_df.groupby(level=0).sum()
                DP_df = DP_df.groupby(level=0).sum()

                # (b) variant_f 压缩到每个 key 一行，顺序对齐 AD_df.index
                #     用 groupby(level=0).first() 保留 ref/alt 等信息时也可扩展
                tmp = variant_f.copy()
                tmp["variant_key"] = variant_key.values
                variant_f_u = tmp.groupby("variant_key", sort=False).first()
                variant_f_u = variant_f_u.loc[AD_df.index].reset_index(drop=True)
                variant_f = variant_f_u

            
            AD_list.append(AD_df)
            DP_list.append(DP_df)

            snv_to_bin_idx = match_snvs_to_bins(variant_f, cnaA.reset_index())
            snv_to_bin_list.append(snv_to_bin_idx)
            print(f"{gt_name } snv_indices min/max:", int(snv_to_bin_idx.min()), int(snv_to_bin_idx.max()))
            print("unique bins in snv_indices:", len(np.unique(snv_to_bin_idx)))

            cnaA = cnaA.loc[:, common_cells]
            cnaB = cnaB.loc[:, common_cells]
            cnaA_list.append(cnaA)
            cnaB_list.append(cnaB)


        for i, name1 in enumerate(tool_names):
            for j, name2 in enumerate(tool_names):
                # if i == j:
                #     continue
                print(f"Calculating LLR under GT={name1} - > Tool1={name1} vs Tool2={name2}")
                print(f"AD1: {AD_list[i].head()}, DP1: {DP_list[i].head()}")
                print(f"AD2: {AD_list[j].head()}, DP2: {DP_list[j].head()}")
                print(f"snv_to_bin1: {snv_to_bin_list[i][:5]}, snv_to_bin2: {snv_to_bin_list[j][:5]}")


                result = calculate_llr(
                    AD_list[i], DP_list[i],
                    AD_list[j], DP_list[j],
                    snv_to_bin_list[i], snv_to_bin_list[j],
                    cnaA_list[i], cnaB_list[i],
                    cnaA_list[j], cnaB_list[j],
                )

                # Attach metadata
                result["GT_Tool"] = name1
                result["Tool1"] = name1
                result["Tool2"] = name2
                result["n_common_cells"] = len(common_cells)
                results.append(result)
                print(result)

        final_df = pd.DataFrame(results)
        out = os.path.join(self.output_dir, "llr_results.csv")
        final_df.to_csv(out, index=False)
        print(f"LLR results saved to {out}")
        return final_df



        AD = scipy.io.mmread(f"{snv_path}/cellSNP.tag.AD.filtered.mtx").tocsr()
        DP = scipy.io.mmread(f"{snv_path}/cellSNP.tag.DP.filtered.mtx").tocsr()

        non_empty_mask = (AD.getnnz(axis=1) > 0) | (DP.getnnz(axis=1) > 0)
        n_total = AD.shape[0]
        n_kept = int(non_empty_mask.sum())
        print(f"Total variants: {n_total}, kept after removing all-zero rows: {n_kept}")

        AD_f = AD[non_empty_mask]
        DP_f = DP[non_empty_mask]
        variant_f = variant_pos_df.loc[non_empty_mask].reset_index(drop=True)
        print("Filtered AD shape:", AD_f.shape)
        print("Filtered DP shape:", DP_f.shape)
        print("Filtered variant_info length:", len(variant_f))

        AD_df = pd.DataFrame.sparse.from_spmatrix(
            AD_f,
            index=variant_f.index,   
            columns=cell_list
        )

        DP_df = pd.DataFrame.sparse.from_spmatrix(
            DP_f,
            index=variant_f.index,   
            columns=cell_list
        )

        results = []

        for path1_A, path1_B, name1 in zip(tool_hap1_files, tool_hap2_files, tool_names):
            for path2_A, path2_B, name2 in zip(tool_hap1_files, tool_hap2_files, tool_names):   
                print(f"Calculating LLR: {name1} vs {name2}")

                cna1_A = pd.read_csv(path1_A,index_col=0)
                cna1_B = pd.read_csv(path1_B,index_col=0)
                snv_to_bin_idx1 = match_snvs_to_bins(variant_f, cna1_A.reset_index())

                cna2_A = pd.read_csv(path2_A,index_col=0)
                cna2_B = pd.read_csv(path2_B,index_col=0)
                snv_to_bin_idx2 = match_snvs_to_bins(variant_f, cna2_A.reset_index())

                result = calculate_llr(AD_df, DP_df,
                            snv_to_bin_idx1, snv_to_bin_idx2,cna1_A, cna1_B,
                            cna2_A, cna2_B)
                
                result['Tool1'] = name1
                result['Tool2'] = name2

                results.append(result)
        
        final_df = pd.DataFrame(results)
        out = os.path.join(self.output_dir, "llr_results.csv")
        final_df.to_csv(out,index=False)
        print(f"LLR results saved to {out}")
        return final_df

                
    def cnclass(
        self,
        tool_hap1_cna_files: List[str],
        tool_hap2_cna_files: List[str],
        tool_names: List[str],
        type: str = "hcCNA",
        outfile: str = "cnclass_results.csv",
    ) -> pd.DataFrame:

        if type not in ["acCNA", "hcCNA"]:
            raise ValueError("type must be 'acCNA' or 'hcCNA'")
        hap_list = ["minor", "major"] if type == "acCNA" else ["hap1", "hap2"]

        conditions = [
            (">=2", "CN_Gain"), ("=1", "CN_Neutral"), ("=0", "CN_Loss"),
            ("=2", "CN_equal_2"), ("=3", "CN_equal_3"), ("=4", "CN_equal_4"),
            ("=5", "CN_equal_5"), ("=6", "CN_equal_6"), ("=7", "CN_equal_7"),
            ("=8", "CN_equal_8"), ("=9", "CN_equal_9"),("=10", "CN_equal_10")
        ]

        all_rows = []
        for g_h1, g_h2, tool1 in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
            gt_h1_raw = read_and_drop_empty(g_h1)
            gt_h2_raw = read_and_drop_empty(g_h2)
            for f_h1, f_h2, tool2 in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
                p_h1 = read_and_drop_empty(f_h1)
                p_h2 = read_and_drop_empty(f_h2)
                gt_h1, p_h1 =  align_cna_bins(gt_h1_raw, p_h1)
                gt_h2, p_h2 =  align_cna_bins(gt_h2_raw, p_h2)

                gt_h1.set_index("region",inplace=True)
                gt_h2.set_index("region",inplace=True)
                p_h1.set_index("region",inplace=True)
                p_h2.set_index("region",inplace=True)


                gt_h1, p_h1 =  align(gt_h1, p_h1)
                gt_h2, p_h2 =  align(gt_h2, p_h2)

                p_h1 = p_h1.reset_index().rename(columns={"index": "region"})
                p_h2 = p_h2.reset_index().rename(columns={"index": "region"})
                gt_h1 = gt_h1.reset_index().rename(columns={"index": "region"})
                gt_h2 = gt_h2.reset_index().rename(columns={"index": "region"})

                for cond, folder in conditions:

                    save_dir = os.path.join(self.output_dir, f"{tool1}_{tool2}", folder)
                    os.makedirs(save_dir, exist_ok=True)
                    categorize_and_save(gt_h1, p_h1, save_dir, tool2, cond, hap_list[0])
                    categorize_and_save(gt_h2, p_h2, save_dir, tool2, cond, hap_list[1])

                    res = process_folder_for_metrics(save_dir, tool2, hap_list)
                    print(res)
                    for htype, metrics in res.items():
                        row = {"Type": folder, "Haplotype": htype, "Tool1": tool1, "Tool2": tool2}
                        row.update(metrics)
                        all_rows.append(row)

        df = pd.DataFrame(all_rows)
        out = os.path.join(self.output_dir, outfile)
        df.to_csv(out, index=False)
        return df

        
    def segmentation(
        self,
        tool_cna_files: List[str],
        tool_names: List[str],
        threshold: float = 0.95,
        outprefix: str = "segmentation_",
    ):
        overlap_results = []
        metric_results = []
        # merged_cna_event_count = None

        # annotated_by_tool = {}
        # for f, name in zip(tool_cna_files, tool_names):
        #     tool_cna_df = pd.read_csv(f, index_col=0)

        #     tool_annotated_df = annotate_segments(
        #         tool_cna_df.T,
        #         detect_cnv_type=True,
        #         threshold=threshold
        #     )

        #     tool_annotated_df["region"] = (
        #         tool_annotated_df["Chrom"].astype(str)
        #         + ":"
        #         + tool_annotated_df["Start"].astype(int).astype(str)
        #         + "-"
        #         + tool_annotated_df["End"].astype(int).astype(str)
        #     )

        #     print(f"{name}_annotated_df: {tool_annotated_df.head()}-----")

        #     annotated_by_tool[name] = tool_annotated_df

            # # 2) 每个工具的 event count 只算一次并 merge
            # tool_counts = get_seg_cna_event_num(tool_annotated_df).copy()
            # tool_counts = tool_counts.rename(
            #     columns={
            #         "bin_level_count": f"{name}_bin_level_count",
            #         "seg_level_count": f"{name}_seg_level_count",
            #     }
            # )

            # if merged_cna_event_count is None:
            #     merged_cna_event_count = tool_counts
            # else:
            #     merged_cna_event_count = merged_cna_event_count.merge(
            #         tool_counts, on=["size", "type"], how="outer"
            #     )

        # if merged_cna_event_count is None:
        #     merged_cna_event_count = pd.DataFrame(columns=["size", "type"])

        n = len(tool_names)
        for i in range(n):
            for j in range(i + 1, n):
                df1 = pd.read_csv(tool_cna_files[i])
                df2 = pd.read_csv(tool_cna_files[j])

                tool1_df, tool2_df = align_cna_bins(df1, df2)

                print(f"tool1 : {tool1_df.head()}")
                print(f"tool2 : {tool2_df.head()}")
                
                cols_except_region = [c for c in tool1_df.columns if c != "region"]
                tool1_df = tool1_df.dropna(subset=cols_except_region, how="all")
                cols_except_region = [c for c in tool2_df.columns if c != "region"]
                tool2_df = tool2_df.dropna(subset=cols_except_region, how="all")

                print(f"tool1 after drop : {tool1_df.head()}")
                print(f"tool2 after drop : {tool2_df.head()}")

                name1 = tool_names[i]
                name2 = tool_names[j]

                tool1_annotated_df = annotate_segments(
                    tool1_df.set_index("region").T, detect_cnv_type=True,
                    threshold=threshold
                )

                tool1_annotated_df["region"] = (
                    tool1_annotated_df["Chrom"].astype(str)
                    + ":"
                    + tool1_annotated_df["Start"].astype(int).astype(str)
                    + "-"
                    + tool1_annotated_df["End"].astype(int).astype(str)
                )

                tool2_annotated_df = annotate_segments(
                    tool2_df.set_index("region").T,detect_cnv_type=True,
                    threshold=threshold
                )

                tool2_annotated_df["region"] = (
                    tool2_annotated_df["Chrom"].astype(str)
                    + ":"
                    + tool2_annotated_df["Start"].astype(int).astype(str)
                    + "-"
                    + tool2_annotated_df["End"].astype(int).astype(str)
                )
                
                print(f"tool1 after align : {tool2_annotated_df.head()}")
                print(f"tool2 after align: {tool2_annotated_df.head()}")
                # cols_except_region = [c for c in tool1.columns if c != "region"]
                # tool1 = tool1.dropna(subset=cols_except_region, how="all")

                # cols_except_region = [c for c in tool2.columns if c != "region"]
                # tool2 = tool2.dropna(subset=cols_except_region, how="all")
                
                

                result_overlap = get_segment_overlap_ratio(tool1_annotated_df, tool2_annotated_df)
                result_metric = get_segment_metric_both(tool1_annotated_df, tool2_annotated_df)

                result_overlap["Tool1"] = name1
                result_overlap["Tool2"] = name2
                result_metric["Tool1"] = name1
                result_metric["Tool2"] = name2

                overlap_results.append(result_overlap)
                metric_results.append(result_metric)

        # # 4) 写出结果
        # num_cols = [c for c in merged_cna_event_count.columns if c.endswith("_count")]
        # if len(num_cols) > 0:
        #     merged_cna_event_count[num_cols] = (
        #         merged_cna_event_count[num_cols].fillna(0).astype(int)
        #     )

        overlap_df = pd.concat(overlap_results, ignore_index=True) if overlap_results else pd.DataFrame()
        metric_df = pd.concat(metric_results, ignore_index=True) if metric_results else pd.DataFrame()

        overlap_df.to_csv(os.path.join(self.output_dir, f"{outprefix}overlap.csv"), index=False)
        metric_df.to_csv(os.path.join(self.output_dir, f"{outprefix}metrics.csv"), index=False)
        # merged_cna_event_count.to_csv(
        #     os.path.join(self.output_dir, f"{outprefix}complex_cna_count.csv"),
        #     index=False
        # )

    # def cloneSizebycellprofile(
    #     self,
    #     tool_cna_files: List[str],
    #     tool_names: List[str],
    #     outfile: str = "clone_size_by_cell_profile.csv",
    # ) -> pd.DataFrame:

    #     ref_file, ref_name = tool_cna_files[0], tool_names[0]
    #     ref_df = pd.read_csv(ref_file, index_col=0)
    #     ref_df = ref_df.fillna("1|1")
    #     ref_clone = get_cell_profile_size(ref_df)  

    #     out = pd.DataFrame(index=ref_clone.index.copy())
    #     out.index.name = ref_clone.index.name
    #     out[f"{ref_name}_pred_size"] = ref_clone["cluster_size"]

    #     for f, name in zip(tool_cna_files[1:], tool_names[1:]):
    #         pred_df = pd.read_csv(f, index_col=0)
    #         pred_df = pred_df.fillna("1|1")
    #         pred_clone = get_cell_profile_size(pred_df)

    #         aligned = pred_clone.reindex(out.index)
    #         out[f"{name}_pred_size"] = aligned["cluster_size"]

    #     # 按 cluster_size 分组求均值
    #     size_col = "cluster_size"
    #     mean_cols = [c for c in out.columns if c.endswith("_pred_size")]
    #     mean_cols = [c for c in mean_cols if pd.api.types.is_numeric_dtype(out[c])]

    #     df_mean = (
    #         out.groupby(size_col, as_index=False)[mean_cols]
    #         .mean()
    #         .rename(columns={size_col: "cluster_size"})
    #     )
    #     df_mean.to_csv(os.path.join(self.output_dir, f"mean_{outfile}"), index=False)

    #     return out


    def get_cell_overlap(
        self,
        tool_cna_files: List[str],
        tool_names: List[str],
        outprefix: str = "cell_overlap",
    ) -> pd.DataFrame:
        
        tool_cells = {}
        for path, name in zip(tool_cna_files, tool_names):
            df = read_and_drop_empty(path)

            tool_cells[name] = set(df.columns[1:])

        result_df = pd.DataFrame(index=tool_names, columns=tool_names, dtype="Int64")
        for t1 in tool_names:
            for t2 in tool_names:
                result_df.loc[t1, t2] = len(tool_cells[t1] & tool_cells[t2])

        common_cells = set.intersection(*(tool_cells[t] for t in tool_names)) if tool_names else set()
        n_common = len(common_cells)

        out = os.path.join(self.output_dir, f"{outprefix}_common_cells_{n_common}.csv")
        result_df.to_csv(out)

        print(f"Cell Overlap matrix saved to {out}")

    
    def cloneSizebycluster(
        self,
        tool_cluster_files: List[str],
        tool_names: List[str],
        outfile: str = "clone_size_by_cluster.csv",
    ) -> pd.DataFrame:

        result_list = []
        for f1 , name1 in zip(tool_cluster_files, tool_names):
            gt_df = pd.read_csv(f1, index_col = 0)

            gt_clone = get_cluster_size(gt_df)

            for f2, name2 in zip(tool_cluster_files, tool_names):
                pred_df = pd.read_csv(f2,index_col = 0)

                pred_clone = get_cluster_size(pred_df)

                aligned = pred_clone.reindex(gt_clone.index)

                gt_clone[f"{name2}_pred_size"] = aligned['cluster_size']
                gt_clone['GT_Tool'] = name1
                # gt_clone[f"{name}_pred_cell"] = aligned.index

            result_list.append(gt_clone)

        result_df= pd.concat(result_list)

        result_df.to_csv(os.path.join(self.output_dir, outfile))

        str_cols = ["cluster_size",'GT_Tool']

        mean_cols = result_df.columns.drop(str_cols) 
        mean_cols = [c for c in mean_cols if pd.api.types.is_numeric_dtype(result_df[c])]

        df_mean = (
            result_df
            .groupby(str_cols, as_index=False)[mean_cols]
            .mean()
            # .rename(columns={size_col: "cluster_size"})
        )

        df_mean.to_csv(os.path.join(self.output_dir, f"mean_{outfile}"), index=False)

    def cloneSizebycellprofile(
        self,
        tool_cna_files: List[str],
        tool_names: List[str],
        outfile: str = "clone_size_by_cell_profile.csv",
    ) -> pd.DataFrame:

        result_list = []
        for f1 , name1 in zip(tool_cna_files, tool_names):
            gt_df = pd.read_csv(f1, index_col = 0)
            gt_df.fillna("1|1",inplace=True)

            gt_clone = get_cell_profile_size(gt_df)

            for f2, name2 in zip(tool_cna_files, tool_names):
                pred_df = pd.read_csv(f2,index_col = 0)
                pred_df.fillna("1|1",inplace=True)

                pred_clone = get_cell_profile_size(pred_df)

                aligned = pred_clone.reindex(gt_clone.index)

                gt_clone[f"{name2}_pred_size"] = aligned['cluster_size']
                gt_clone['GT_Tool'] = name1
                # gt_clone[f"{name}_pred_cell"] = aligned.index

            result_list.append(gt_clone)

        result_df= pd.concat(result_list)

        result_df.to_csv(os.path.join(self.output_dir, outfile))

        str_cols = ["cluster_size",'GT_Tool']

        mean_cols = result_df.columns.drop(str_cols) 
        mean_cols = [c for c in mean_cols if pd.api.types.is_numeric_dtype(result_df[c])]

        df_mean = (
            result_df
            .groupby(str_cols, as_index=False)[mean_cols]
            .mean()
            # .rename(columns={size_col: "cluster_size"})
        )

        df_mean.to_csv(os.path.join(self.output_dir, f"mean_{outfile}"), index=False)

    def hcPhasing(
        self,
        tool_hap1_cna_files: List[str],
        tool_hap2_cna_files: List[str],
        tool_names: List[str],
        outprefix = "hcPhasing",
        mode = "heterozygous-only",
    ) -> pd.DataFrame:

        rows = []
        for g_h1, g_h2, name1 in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
            g1 = read_and_drop_empty(g_h1)
            g2 = read_and_drop_empty(g_h2)
            
            for f_h1, f_h2, name in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
                t1 = read_and_drop_empty(f_h1)
                t2 = read_and_drop_empty(f_h2)
                print(f"shape: {g1.shape}, hap1 shape: {t1.shape},hap2 shape: {t2.shape}")


                g1_df, t1_df = align_cna_bins(g1.copy(), t1.copy())
                g2_df, t2_df = align_cna_bins(g2.copy(), t2.copy())

                g1_df.set_index("region",inplace=True)
                t1_df.set_index("region",inplace=True)
                g2_df.set_index("region",inplace=True)
                t2_df.set_index("region",inplace=True)

                print(f"After align change shape: {g1_df.shape}, hap1 shape: {t1_df.shape},hap2 shape: {t2_df.shape}")

                g1_bin, _ = phase_to_binary(g1_df, g2_df)

                t1_bin, _ = phase_to_binary(t1_df, t2_df)

                mismatch_error,_ = real_cell_mismatch_error(t1_bin, g1_bin, mode)
                switch_error,_ = real_cell_switch_error(t1_bin, g1_bin, mode)

                result = pd.DataFrame({
                    "mismatch_ratio": [mismatch_error],
                    "switch_error_ratio": [switch_error]
                })

                result["tool_name"] = name

                result['GT_tool'] = name1
                rows.extend(result.to_dict("records"))

        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(self.output_dir, f"{outprefix}_{mode}.csv"), index=False)

        return df




        
