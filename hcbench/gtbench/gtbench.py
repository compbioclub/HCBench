# hcbench/gtbench/gtbench.py
from math import nan
from operator import gt
import os
import ast
import subprocess
import numpy as np
import pandas as pd
from typing import List, Tuple, Dict, Optional

from scipy.stats import spearmanr
from sklearn.metrics import (
    mean_squared_error, roc_auc_score, precision_recall_curve, auc,
    accuracy_score, precision_score, recall_score, f1_score,
    confusion_matrix, cohen_kappa_score, brier_score_loss
)
from Bio import Phylo

from hcbench.utils import align, align_cna_bins, annotate_segments, categorize_and_save, compute_rd_cn_l1, evaluate_NA_ratio, evaluate_haplotype_predictions, fast_mode_any, find_mirrored_clones, get_cell_profile_size, get_cluster_size, get_seg_cna_event_num, get_seg_mirror_subclonal, get_segment_metric, get_segment_overlap_ratio, process_folder_for_metrics, process_folder_for_metrics_clone,fast_mode_num, read_and_drop_empty
from sklearn.metrics import adjusted_rand_score as adjustedRandIndex, mean_squared_error
from sklearn.metrics import adjusted_mutual_info_score as AMI
from hcbench.parsers.utils import split_all_regions
from ..utils import count_unique_cells
from ..realbench.utils import create_lazac_input,get_final_parsimony_score

class GTBench:

    def __init__(self, output_dir: str = "./output"):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)


    @staticmethod
    def _normalize_hap_label(s: pd.Series) -> pd.Series:

        mapping = {
            'maternal': 'hap1', 'paternal': 'hap2',
            'Maternal': 'hap1', 'Paternal': 'hap2',
            'MATERNAL': 'hap1', 'PATERNAL': 'hap2',
            'hap1': 'hap1', 'hap2': 'hap2', 'HAP1': 'hap1', 'HAP2': 'hap2'
        }
        return s.map(mapping).fillna(s)

    @staticmethod
    def _is_hap(series_val: str, target: str) -> bool:
        val = str(series_val)
        if val.lower() == 'maternal':
            val = 'hap1'
        elif val.lower() == 'paternal':
            val = 'hap2'
        return val == target

    # ========= 1) bin-level: cndetect =========
    def cndetect(
        self,
        tool_cna_files: List[str],
        cna_profile_file: str,
        tool_names: List[str],
        haplotype: str = "combined",
        profile_bin_size = 100000,
        outfile: str = "bin_level_results.csv",
        index_col: Optional[int] = 0,
    ) -> pd.DataFrame:

        truth_df = read_and_drop_empty(cna_profile_file)
        truth_df.set_index("region",inplace=True)
        results = []

        for path, name in zip(tool_cna_files,tool_names):
            print(name)
            pred = read_and_drop_empty(path)

            print(f"gt shape: {truth_df.shape}, {name} shape: {pred.shape}")
            pred = split_all_regions(pred.set_index("region"), profile_bin_size)
            # pred = pred.reset_index().rename(columns={"index": "region"})

            truth,pred = align(truth_df,pred)
            print(f"after align gt shape: {truth.shape}, {name} shape: {pred.shape}")

            rmse, scc, acc = evaluate_haplotype_predictions(pred, truth, haplotype)
            results.append({"Tool": name, "RMSE": rmse, "SCC": scc, "ACC": acc})
            print(results)

        df = pd.DataFrame(results)
        out = os.path.join(self.output_dir, outfile)
        df.to_csv(out, index=False)
        return df


    def cnclass(
        self,
        tool_hap1_cna_files: List[str],
        tool_hap2_cna_files: List[str],
        tool_names: List[str],
        profile_hap1_cna_file: str,
        profile_hap2_cna_file: str,
        type: str = "hcCNA",
        profile_bin_size = 100000,
        outfile: str = "cnclass_results.csv",
    ) -> pd.DataFrame:

        if type not in ["acCNA", "hcCNA"]:
            raise ValueError("type must be 'acCNA' or 'hcCNA'")
        hap_list = ["minor", "major"] if type == "acCNA" else ["hap1", "hap2"]

        gt_h1_r = read_and_drop_empty(profile_hap1_cna_file)
        gt_h1_r.set_index("region",inplace=True)
        gt_h2_r = read_and_drop_empty(profile_hap2_cna_file)
        gt_h2_r.set_index("region",inplace=True)


        conditions = [
            (">=2", "CN_Gain"), ("=1", "CN_Neutral"), ("=0", "CN_Loss"),
            ("=2", "CN_equal_2"), ("=3", "CN_equal_3"), ("=4", "CN_equal_4"),
            ("=5", "CN_equal_5"), ("=6", "CN_equal_6"), ("=7", "CN_equal_7"),
            ("=8", "CN_equal_8"), ("=9", "CN_equal_9"),("=10", "CN_equal_10")
        ]

        all_rows = []
        for f_h1, f_h2, tool in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
            p_h1 = pd.read_csv(f_h1).fillna(-1) 
            p_h2 = pd.read_csv(f_h2).fillna(-1)

            p_h1 = split_all_regions(p_h1.set_index("region"), profile_bin_size)
            p_h2 = split_all_regions(p_h2.set_index("region"), profile_bin_size)

            gt_h1, p_h1 =  align(gt_h1_r, p_h1)
            gt_h2, p_h2 =  align(gt_h2_r, p_h2)

            p_h1 = p_h1.reset_index().rename(columns={"index": "region"})
            p_h2 = p_h2.reset_index().rename(columns={"index": "region"})

            gt_h1 = gt_h1.reset_index().rename(columns={"index": "region"})
            gt_h2 = gt_h2.reset_index().rename(columns={"index": "region"})


            for cond, folder in conditions:

                save_dir = os.path.join(self.output_dir, tool, folder)
                os.makedirs(save_dir, exist_ok=True)

                categorize_and_save(gt_h1, p_h1, save_dir, tool, cond, hap_list[0])
                categorize_and_save(gt_h2, p_h2, save_dir, tool, cond, hap_list[1])

                res = process_folder_for_metrics_clone(save_dir, tool, hap_list)
                for htype, clones in res.items():
                    for clone, metrics in clones.items():
                        row = {"Type": folder, "Haplotype": htype, "Clone": clone, "Tool": tool}
                        row.update(metrics)
                        all_rows.append(row)

        df = pd.DataFrame(all_rows)
        out = os.path.join(self.output_dir, outfile)
        df.to_csv(out, index=False)
        return df

    @staticmethod
    def _align_for_csv(gdf: pd.DataFrame, pdf: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        cols = gdf.columns.intersection(pdf.columns)

        return gdf.loc[:, cols], pdf.loc[:, cols]

    # ========= 3) hccnchange =========
    def hccnchange(
        self,
        tool_hap1_cna_files: List[str],
        tool_hap2_cna_files: List[str],
        tool_names: List[str],
        changes_file: str,
        profile_bin_size = 100000,
        outfile: str = "evolution_onset_CN_Change.csv",
    ) -> pd.DataFrame:

        change_truth = read_and_drop_empty(changes_file)
    
        if 'Haplotype' in change_truth.columns:
            change_truth['Haplotype'] = self._normalize_hap_label(change_truth['Haplotype'])

        change_truth = split_all_regions(change_truth.set_index("Segment"), profile_bin_size)
        change_truth = change_truth.reset_index().rename(columns={"index": "Segment"})

        results = []

        for f_h1, f_h2, name in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
            p_h1 = read_and_drop_empty(f_h1)
            p_h2 = read_and_drop_empty(f_h2)

            p_h1 = split_all_regions(p_h1.set_index("region"), profile_bin_size)
            p_h1 = p_h1.reset_index().rename(columns={"index": "region"})
            p_h2 = split_all_regions(p_h2.set_index("region"), profile_bin_size)
            p_h2 = p_h2.reset_index().rename(columns={"index": "region"})


            h1 = self._onset_join(change_truth, p_h1, 'hap1')
            h2 = self._onset_join(change_truth, p_h2, 'hap2')

            comb = pd.concat([h1, h2], ignore_index=True)

            for t in comb['Type'].unique().tolist():

                gt = comb[comb['Type'] == t]['Change']
                pd_ = comb[comb['Type'] == t]['Change_predict']
                rmse = self._cn_change_rmse(gt, pd_)
                acc  = self._cn_change_acc(gt, pd_)
                results.append({"Tool": name, "Type": t, "RMSE": rmse, "ACC": acc})

        df = pd.DataFrame(results)
        df.to_csv(os.path.join(self.output_dir, outfile), index=False)
        return df


    @staticmethod
    def _cn_change_acc(gt: pd.Series, pd_: pd.Series) -> float:
        p  = gt.str.split('->').str[0].apply(pd.to_numeric, errors='coerce')
        pp = pd_.str.split('->').str[0].apply(pd.to_numeric, errors='coerce')
        c  = gt.str.split('->').str[1].apply(pd.to_numeric, errors='coerce')
        cp = pd_.str.split('->').str[1].apply(pd.to_numeric, errors='coerce')
        mask = ~p.isna() & ~pp.isna() & ~c.isna() & ~cp.isna()
        if not mask.any():
            return float("nan")
        return float(((p - c) == (pp - cp))[mask].mean())

    @staticmethod
    def _cn_change_rmse(gt: pd.Series, pd_: pd.Series) -> float:
        p  = gt.str.split('->').str[0].apply(pd.to_numeric, errors='coerce')
        pp = pd_.str.split('->').str[0].apply(pd.to_numeric, errors='coerce')
        c  = gt.str.split('->').str[1].apply(pd.to_numeric, errors='coerce')
        cp = pd_.str.split('->').str[1].apply(pd.to_numeric, errors='coerce')
        diff_true = (p - c).astype(float)
        diff_pred = (pp - cp).astype(float)
        mask = ~diff_true.isna() & ~diff_pred.isna()
        if not mask.any():
            return float("nan")
        return float(np.sqrt(mean_squared_error(diff_true[mask], diff_pred[mask])))

    # ========= 4) hccnstable =========
    def hccnstable(
        self,
        tool_hap1_cna_files: List[str],
        tool_hap2_cna_files: List[str],
        tool_names: List[str],
        changes_file: str,
        tree_file: str,
        outfile: str = "evolution_cn_stability_acc.csv",
        profile_bin_size = 100000
    ) -> pd.DataFrame:

        tree = Phylo.read(tree_file, "newick")
        change_df = read_and_drop_empty(changes_file)

        change_df = self._add_check_list(tree, change_df)

        change_df = split_all_regions(change_df.set_index("Segment"), profile_bin_size)
        change_df = change_df.reset_index().rename(columns={"index": "Segment"})

        if 'Haplotype' in change_df.columns:
            change_df['Haplotype'] = self._normalize_hap_label(change_df['Haplotype'])
        

        results = []
        for f_h1, f_h2, name in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
            p_h1 = read_and_drop_empty(f_h1)
            p_h2 = read_and_drop_empty(f_h2)

            p_h1 = split_all_regions(p_h1.set_index("region"), profile_bin_size)
            p_h1 = p_h1.reset_index().rename(columns={"index": "region"})
            p_h2 = split_all_regions(p_h2.set_index("region"), profile_bin_size)
            p_h2 = p_h2.reset_index().rename(columns={"index": "region"})

            proc = change_df.copy()

            self._process_change_rows(proc, p_h1, p_h2)

            for t in proc['Type'].unique().tolist():
                acc = proc[proc['Type'] == t]['result'].mean()
                results.append({"Tool": name, "Type": t, "ACC": float(acc)})

        df = pd.DataFrame(results)
        df.to_csv(os.path.join(self.output_dir, outfile), index=False)
        return df


    @staticmethod
    def _is_desc_or_self(tree, anc, des) -> bool:
        ancestor = tree.find_any(name=anc); descendant = tree.find_any(name=des)
        if ancestor is None or descendant is None: return False
        if ancestor == descendant: return True
        for cl in ancestor.find_clades():
            if cl == descendant:
                return True
        return False

    @staticmethod
    def _get_desc(tree, anc) -> List[str]:
        node = tree.find_any(name=anc)
        if node is None:
            return []
        return [cl.name for cl in node.find_clades() if cl.name is not None]

    def _add_check_list(self, tree, change_df: pd.DataFrame) -> pd.DataFrame:
        if 'check_child_list' not in change_df.columns:
            change_df['check_child_list'] = None
        change_df = change_df.astype({'check_child_list': 'object'})
        for idx, row in change_df.iterrows():
            start = row['Child']; segment = row['Segment']; hap = row['Haplotype']
            descendant = self._get_desc(tree, start)
            sub_clone = change_df[(change_df['Segment'] == segment) & (change_df['Haplotype'] == hap)]['Child'].unique()
            dels = []
            for cl in sub_clone:
                if cl == start:
                    continue
                elif self._is_desc_or_self(tree, start, cl):
                    dels += list(self._get_desc(tree, cl))
            update = list(set(descendant) - set(dels))
            change_df.at[idx, 'check_child_list'] = update
        return change_df

    @staticmethod
    def _mode_for_prefix_row(row: pd.Series, prefix: str, df_row: pd.Series, full_df: pd.DataFrame):
        cols = [c for c in full_df.columns if c.startswith(f"{prefix}_")]
        if not cols:
            return None
        m = df_row[cols].mode()
        return m.iloc[0] if not m.empty else None

    def _process_change_rows(self, change_df: pd.DataFrame, p_h1: pd.DataFrame, p_h2: pd.DataFrame):
        change_df['result'] = 0
        for idx, row in change_df.iterrows():
            region = row['Segment']
            hap = row['Haplotype']
            target = p_h1 if self._is_hap(hap, 'hap1') else p_h2
            tool_row = target[target['region'] == region]
            if tool_row.empty:
                continue
            tool_row = tool_row.iloc[0]
            clone_list = row['check_child_list']
            if isinstance(clone_list, str):
                clone_list = ast.literal_eval(clone_list)
            ok = True
            for cl in (clone_list or []):
                v = self._mode_for_prefix_row(row, cl, tool_row, target)
                if v != int(row['Change'].split('->')[1]):
                    ok = False; break
            change_df.at[idx, 'result'] = 1 if ok else 0

    # ========= 5) hconsetacc =========
    def hconsetacc(
        self,
        tool_hap1_cna_files: List[str],
        tool_hap2_cna_files: List[str],
        tool_names: List[str],
        changes_file: str,
        profile_bin_size = 100000,
        outfile: str = "evolution_onset_acc.csv",
    ) -> pd.DataFrame:

        change_truth = read_and_drop_empty(changes_file)

        change_truth = split_all_regions(change_truth.set_index("Segment"), profile_bin_size)
        change_truth = change_truth.reset_index().rename(columns={"index": "Segment"})


        if 'Haplotype' in change_truth.columns:
            change_truth['Haplotype'] = self._normalize_hap_label(change_truth['Haplotype'])
        results = []

        for f_h1, f_h2, name in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):            
            p_h1 = read_and_drop_empty(f_h1); p_h2 = read_and_drop_empty(f_h2)

            p_h1 = split_all_regions(p_h1.set_index("region"), profile_bin_size)
            p_h1 = p_h1.reset_index().rename(columns={"index": "region"})
            p_h2 = split_all_regions(p_h2.set_index("region"), profile_bin_size)
            p_h2 = p_h2.reset_index().rename(columns={"index": "region"})


            h1 = self._onset_join(change_truth, p_h1, 'hap1')
            h2 = self._onset_join(change_truth, p_h2, 'hap2')
            comb = pd.concat([h1, h2], ignore_index=True)

            comb.to_csv(os.path.join(self.output_dir, f"{name}_comb_combined.csv"))


            for t in comb['Type'].unique().tolist():
                gt = comb[comb['Type'] == t]['Change']
                pd_ = comb[comb['Type'] == t]['Change_predict']

                mask = gt.notna() & pd_.notna()
                acc = float((pd_[mask] == gt[mask]).mean()) if mask.any() else float("nan")
                results.append({"Tool": name, "Type": t, "ACC": acc})

        df = pd.DataFrame(results)
        df.to_csv(os.path.join(self.output_dir, outfile), index=False)
        return df

   
    @staticmethod
    def _onset_join(change_df: pd.DataFrame, pred_df: pd.DataFrame, hap: str) -> pd.DataFrame:
        filt = change_df[GTBench._normalize_hap_label(change_df['Haplotype']) == hap]
        merged = pred_df.merge(filt, left_on='region', right_on='Segment', how='inner')
        if merged.empty:
            return pd.DataFrame(columns=['Parent','Child','Haplotype','Type','Segment','Change','Parent_predict_num','Child_predict_num','Change_predict'])

        parent_vals = merged["Parent"].astype(str).to_numpy()
        child_vals  = merged["Child"].astype(str).to_numpy()

        all_cols = merged.columns.tolist()

        def build_prefix_map(prefixes):
            uniq = np.unique(prefixes)
            return {
                p: np.array([i for i, c in enumerate(all_cols) if c.startswith(p)], dtype=int)
                for p in uniq
            }

        parent_colmap = build_prefix_map(parent_vals)
        child_colmap  = build_prefix_map(child_vals)

        arr = merged.to_numpy()      # shape = (n_rows, n_cols)
        n = arr.shape[0]

        parent_pred = np.empty(n, dtype=object)
        child_pred  = np.empty(n, dtype=object)

        for i in range(n):
            p = parent_vals[i]
            cols = parent_colmap.get(p)
            if cols is None or len(cols) == 0:
                parent_pred[i] = None
            else:
                parent_pred[i] = fast_mode_num(arr[i, cols].astype(float))

            c = child_vals[i]
            cols = child_colmap.get(c)
            if cols is None or len(cols) == 0:
                child_pred[i] = None
            else:
                child_pred[i] = fast_mode_num(arr[i, cols].astype(float))

        out = merged[[
            'Parent','Child','Haplotype','Type','Segment','Change'
        ]].copy()

        out['Parent_predict_num'] = parent_pred
        out['Child_predict_num']  = child_pred

        m = out["Parent_predict_num"].notna() & out["Child_predict_num"].notna()

        out["Change_predict"] = pd.NA
        out.loc[m, "Change_predict"] = (
            out.loc[m, "Parent_predict_num"].astype(int).astype(str)
            + "->" +
            out.loc[m, "Child_predict_num"].astype(int).astype(str)
        )
        return out


    # ========= 6) hconsetcn =========
    def hconsetcn(
        self,
        tool_hap1_cna_files: List[str],
        tool_hap2_cna_files: List[str],
        tool_names: List[str],
        changes_file: str,
        outfile: str = "evolution_onset_parent_CN.csv",
        profile_bin_size = 100000
    ) -> pd.DataFrame:

        change_truth_r = read_and_drop_empty(changes_file)
        if 'Haplotype' in change_truth_r.columns:
            change_truth_r['Haplotype'] = self._normalize_hap_label(change_truth_r['Haplotype'])
        results = []

        change_truth_r = split_all_regions(change_truth_r.set_index("Segment"), profile_bin_size)
        change_truth_r = change_truth_r.reset_index().rename(columns={"index": "Segment"})

        for f_h1, f_h2, name in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
            p_h1 = read_and_drop_empty(f_h1); p_h2 = read_and_drop_empty(f_h2)

            p_h1 = split_all_regions(p_h1.set_index("region"), profile_bin_size)
            p_h1 = p_h1.reset_index().rename(columns={"index": "region"})
            p_h2 = split_all_regions(p_h2.set_index("region"), profile_bin_size)
            p_h2 = p_h2.reset_index().rename(columns={"index": "region"})

            for t in change_truth_r['Type'].unique().tolist():
                change_truth = change_truth_r[change_truth_r['Type'] == t]
    
                h1 = self._onset_join(change_truth, p_h1, 'hap1')
                h2 = self._onset_join(change_truth, p_h2, 'hap2')
                comb = pd.concat([h1, h2], ignore_index=True)

                print(f"{change_truth['Type'].unique().tolist()}: {change_truth['Type'].shape}")

                gt = comb['Change']
                pd_p = comb["Parent_predict_num"]
                pd_c = comb["Child_predict_num"]
                acc = self._parent_onset_acc(gt, pd_p)
                pr, _ = self._parent_child_rmse(gt, pd_p, pd_c)
                results.append({"Tool": name, "Type": t, "RMSE": pr, "ACC": acc})

        df = pd.DataFrame(results)
        df.to_csv(os.path.join(self.output_dir, outfile), index=False)
        return df

    @staticmethod
    def _parent_onset_acc(gt: pd.Series, pd_: pd.Series) -> float:
        p  = gt.str.split('->').str[0].apply(pd.to_numeric, errors='coerce')
        pp = pd.to_numeric(pd_, errors="coerce")
        mask = ~p.isna() & ~pp.isna()
        if not mask.any():
            return float("nan")
        return float((p[mask] == pp[mask]).mean())

    @staticmethod
    def _parent_child_rmse(gt: pd.Series, pd_p: pd.Series, pd_c: pd.Series) -> Tuple[float, float]:
        p  = gt.str.split('->').str[0].apply(pd.to_numeric, errors='coerce')
        pp = pd.to_numeric(pd_p, errors="coerce")
        c  = gt.str.split('->').str[1].apply(pd.to_numeric, errors='coerce')
        cp = pd.to_numeric(pd_c, errors="coerce")
        m1 = ~p.isna() & ~pp.isna()
        m2 = ~c.isna() & ~cp.isna()
        pr = float(np.sqrt(mean_squared_error(p[m1], pp[m1]))) if m1.any() else float("nan")
        cr = float(np.sqrt(mean_squared_error(c[m2], cp[m2]))) if m2.any() else float("nan")
        return pr, cr

    # ========= 7) hcPhasing =========
    # def hcPhasing(
    #     self,
    #     tool_hap1_cna_files: List[str],
    #     tool_hap2_cna_files: List[str],
    #     tool_names: List[str],
    #     ground_truth_hap1_file: str,
    #     ground_truth_hap2_file: str,
    #     outprefix = "hcPhasing",
    #     profile_bin_size = 100000,
    #     mask_both = True,
    #     output_all = False
    # ) -> pd.DataFrame:

    #     print("read gt")
    #     g1_r = read_and_drop_empty(ground_truth_hap1_file)
    #     g2_r = read_and_drop_empty(ground_truth_hap2_file)

    #     g1_r.set_index("region",inplace=True)
    #     g2_r.set_index("region",inplace=True)

    #     rows = []
    #     for f_h1, f_h2, name in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
    #         t1 = read_and_drop_empty(f_h1)
    #         t2 = read_and_drop_empty(f_h2)

    #         t1 = split_all_regions(t1.set_index("region"), profile_bin_size)
    #         t2 = split_all_regions(t2.set_index("region"), profile_bin_size)
    #         g1, t1 =  align(g1_r, t1)
    #         g2, t2 =  align(g2_r, t2)

    #         print(f"After align change shape: {g1.shape}, hap1 shape: {t1.shape},hap2 shape: {t2.shape}")

    #         g1_bin, g2_bin = self._phase_to_binary(g1, g2)

    #         h1_bin, h2_bin = self._phase_to_binary(t1, t2)

    #         eval_fn = self._eval_mismatch_switch_both if mask_both else self._eval_mismatch_switch_gt

    #         part_a = eval_fn(g1_bin, g2_bin, h1_bin, h2_bin, name)
    #         part_b = eval_fn(g2_bin, g1_bin, h1_bin, h2_bin, name)

    #         df_a = pd.DataFrame(part_a)
    #         df_b = pd.DataFrame(part_b)

    #         if df_a.empty and df_b.empty:
    #             continue
    #         if df_a.empty:
    #             rows.extend(part_b)
    #             continue
    #         if df_b.empty:
    #             rows.extend(part_a)
    #             continue

    #         for c in ["mismatch_count", "mismatch_ratio", "total",
    #               "switch_error_count", "total_switch_compare_count", "switch_error_ratio"]:
    #             if c in df_a.columns:
    #                 df_a[c] = pd.to_numeric(df_a[c], errors="coerce")
    #             if c in df_b.columns:
    #                 df_b[c] = pd.to_numeric(df_b[c], errors="coerce")

    #         total_a = float(df_a["mismatch_count"].sum())
    #         total_b = float(df_b["mismatch_count"].sum())

    #         best = df_b if total_b < total_a else df_a
    #         rows.extend(best.to_dict("records"))

    #     df = pd.DataFrame(rows)
    #     df.to_csv(os.path.join(self.output_dir, f"{outprefix}_perclone.csv"), index=False)

    #     if output_all and not df.empty:
    #         num_cols = ["mismatch_count", "total", "switch_error_count", "total_switch_compare_count"]
    #         for c in num_cols:
    #             if c in df.columns:
    #                 df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    #         tool_sum = (
    #             df.groupby("tool_name", as_index=False)[
    #                 ["mismatch_count", "total", "switch_error_count", "total_switch_compare_count"]
    #             ].sum()
    #         )

    #         tool_sum["mismatch_ratio"] = np.where(
    #             tool_sum["total"] > 0,
    #             tool_sum["mismatch_count"] / tool_sum["total"],
    #             np.nan
    #         )
    #         tool_sum["switch_error_ratio"] = np.where(
    #             tool_sum["total_switch_compare_count"] > 0,
    #             tool_sum["switch_error_count"] / tool_sum["total_switch_compare_count"],
    #             np.nan
    #         )

    #         tool_sum = tool_sum[
    #             ["tool_name",
    #             "mismatch_count", "total", "mismatch_ratio",
    #             "switch_error_count", "total_switch_compare_count", "switch_error_ratio"]
    #         ]
    #         tool_sum.to_csv(os.path.join(self.output_dir, f"{outprefix}_total.csv"), index=False)

    #     return df

    # ========= 7) hcPhasing =========
    def hcPhasing(
        self,
        tool_hap1_cna_files: List[str],
        tool_hap2_cna_files: List[str],
        tool_names: List[str],
        ground_truth_hap1_file: str,
        ground_truth_hap2_file: str,
        outprefix = "hcPhasing",
        profile_bin_size = 100000,
        mask_both = True,
        output_all = False
    ) -> pd.DataFrame:

        print("read gt")
        g1_r = read_and_drop_empty(ground_truth_hap1_file)
        g2_r = read_and_drop_empty(ground_truth_hap2_file)

        g1_r.set_index("region",inplace=True)
        g2_r.set_index("region",inplace=True)

        rows = []
        for f_h1, f_h2, name in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
            t1 = read_and_drop_empty(f_h1)
            t2 = read_and_drop_empty(f_h2)

            t1 = split_all_regions(t1.set_index("region"), profile_bin_size)
            t2 = split_all_regions(t2.set_index("region"), profile_bin_size)
            g1, t1 =  align(g1_r, t1)
            g2, t2 =  align(g2_r, t2)

            print(f"After align change shape: {g1.shape}, hap1 shape: {t1.shape},hap2 shape: {t2.shape}")

            g1_bin, g2_bin = self._phase_to_binary(g1, g2)

            h1_bin, h2_bin = self._phase_to_binary(t1, t2)

            eval_fn = self._eval_mismatch_switch_both if mask_both else self._eval_mismatch_switch_gt

            
            # tool_sum.to_csv(os.path.join(self.output_dir, f"{outprefix}_total.csv"), index=False)

        return df

    @staticmethod
    def _phase_to_binary(hap1_df: pd.DataFrame, hap2_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
       
        h1 = (hap1_df > hap2_df).astype(int)
        h1 = h1.mask(hap1_df == hap2_df, -1).mask((hap1_df.isna() | hap2_df.isna()), -1)
        h2 = (hap2_df > hap1_df).astype(int)
        h2 = h2.mask(hap1_df == hap2_df, -1).mask((hap1_df.isna() | hap2_df.isna()), -1)
        return h1, h2

    @staticmethod
    def _eval_mismatch_switch_gt(g1: pd.DataFrame, g2: pd.DataFrame,
                              h1: pd.DataFrame, h2: pd.DataFrame, tool_name: str) -> List[Dict]:
        # g1 = g1.dropna(axis=1, how='all'); g2 = g2.dropna(axis=1, how='all')
   
        # idx = g1.index.intersection(h1.index)
        # col = g1.columns.intersection(h1.columns)
        # g1, g2, h1, h2 = g1.loc[idx, col], g2.loc[idx, col], h1.loc[idx, col], h2.loc[idx, col]

        col = g1.columns

        prefixes = col.str.extract(r"(^[^_]+)_")[0].unique()
        results = []

        for pref in prefixes:
            grp_cols = [c for c in col if c.startswith(f"{pref}_")]
            g1_g, g2_g = g1[grp_cols], g2[grp_cols]
            h1_g,  h2_g = h1[grp_cols],  h2[grp_cols]

            # mismatch
            # mask_1 = (g1_g != -1) & (h1_g != -1)
            # mask_2 = (g2_g != -1) & (h2_g != -1)

            mask_1 = g1_g != -1; mask_2 = g2_g != -1

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

            valid1 = (g1_cur != -1) & (g1_next != -1)
            # hap2 相邻 bin 同时有效
            valid2 = (g2_cur != -1) & (g2_next != -1)


            sw1_mask = ((g1_cur != h1_cur) | (g1_next != h1_next)) & valid1
            sw2_mask = ((g2_cur != h2_cur) | (g2_next != h2_next)) & valid2

            total_sw = valid1.sum()
            sw_1 = sw1_mask.sum()

            results.append({
                "tool_name": tool_name,
                "cell_group": pref,
                "mismatch_count": int(mm),
                "total": int(total_m),
                "mismatch_ratio": (mm / total_m) if total_m > 0 else None,
                "switch_error_count": int(sw_1),
                "total_switch_compare_count": int(total_sw),
                "switch_error_ratio": (sw_1 / total_sw) if total_sw > 0 else None,
            })
        return results
    
    @staticmethod
    def _eval_mismatch_switch_both(g1: pd.DataFrame, g2: pd.DataFrame,
                              h1: pd.DataFrame, h2: pd.DataFrame, tool_name: str) -> List[Dict]:
        # g1 = g1.dropna(axis=1, how='all'); g2 = g2.dropna(axis=1, how='all')
   
        # idx = g1.index.intersection(h1.index)
        # col = g1.columns.intersection(h1.columns)
        # g1, g2, h1, h2 = g1.loc[idx, col], g2.loc[idx, col], h1.loc[idx, col], h2.loc[idx, col]

        col = g1.columns

        prefixes = col.str.extract(r"(^[^_]+)_")[0].unique()
        results = []

        for pref in prefixes:
            grp_cols = [c for c in col if c.startswith(f"{pref}_")]
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
                "cell_group": pref,
                "mismatch_count": int(mm),
                "total": int(total_m),
                "mismatch_ratio": (mm / total_m) if total_m > 0 else None,
                "switch_error_count": int(sw_1),
                "total_switch_compare_count": int(total_sw),
                "switch_error_ratio": (sw_1 / total_sw) if total_sw > 0 else None,
            })
        return results

    # ========= 8) mirrorsubclone =========
    def mirrorsubclone(
        self,
        tool_cna_files: List[str],
        tool_names: List[str],
        changes_file: str,
        profile_bin_size =100000,
        outfile: str = "mirror_subclone_result.csv",
    ) -> pd.DataFrame:

        change_df_r = read_and_drop_empty(changes_file)

        change_df_r['region'] = (
            change_df_r['Chromosome'].astype(str) + ":" +
            change_df_r['Start'].astype(str) + "-" +
            change_df_r['End'].astype(str)
        )

        change_truth_r = split_all_regions(change_df_r.set_index("region"), profile_bin_size)
        change_truth_r = change_truth_r.reset_index().rename(columns={"index": "region"})

        rows = []
        for path, name in zip(tool_cna_files, tool_names):
            pred = read_and_drop_empty(path)

            pred = split_all_regions(pred.set_index("region"), profile_bin_size)
            pred = pred.reset_index().rename(columns={"index": "region"})

            change_df = change_truth_r

            change_df = change_df.dropna()
            combined = self._mirror_merge(change_df, pred)
            print(combined.head())
            if combined.empty:
                rows.append({"Tool": name, "RMSE": float("nan"), "ACC": float("nan")})
                continue

         
            for side in ["Clone1", "Clone2"]:
                combined[f"{side}_hap1_CNA"] = combined[f"{side}_CNA"].astype(str).str.split('|').str[0]
                combined[f"{side}_hap2_CNA"] = combined[f"{side}_CNA"].astype(str).str.split('|').str[1]
                combined[f"{side}_predict_hap1_CNA"] = combined[f"{side}_predict_CNA"].astype(str).str.split('|').str[0]
                combined[f"{side}_predict_hap2_CNA"] = combined[f"{side}_predict_CNA"].astype(str).str.split('|').str[1]


            acc_mask = (
                (combined['Clone1_predict_hap1_CNA'].astype(str) == combined['Clone1_hap1_CNA'].astype(str)) &
                (combined['Clone1_predict_hap2_CNA'].astype(str) == combined['Clone1_hap2_CNA'].astype(str)) &
                (combined['Clone2_predict_hap1_CNA'].astype(str) == combined['Clone2_hap1_CNA'].astype(str)) &
                (combined['Clone2_predict_hap2_CNA'].astype(str) == combined['Clone2_hap2_CNA'].astype(str))
            )
            acc = float(acc_mask.mean())


            def se(a, b):  
                aa = pd.to_numeric(a, errors='coerce').astype(float)
                bb = pd.to_numeric(b, errors='coerce').astype(float)
                return (aa - bb) ** 2

            combined['Clone1_hap1_error'] = se(combined['Clone1_predict_hap1_CNA'], combined['Clone1_hap1_CNA'])
            combined['Clone1_hap2_error'] = se(combined['Clone1_predict_hap2_CNA'], combined['Clone1_hap2_CNA'])
            combined['Clone2_hap1_error'] = se(combined['Clone2_predict_hap1_CNA'], combined['Clone2_hap1_CNA'])
            combined['Clone2_hap2_error'] = se(combined['Clone2_predict_hap2_CNA'], combined['Clone2_hap2_CNA'])

            combined['RMSE_result'] = (
                np.sqrt(combined['Clone1_hap1_error']) +
                np.sqrt(combined['Clone1_hap2_error']) +
                np.sqrt(combined['Clone2_hap1_error']) +
                np.sqrt(combined['Clone2_hap2_error'])
            )
            rmse = float(combined['RMSE_result'].mean())

            rows.append({"Tool": name, "RMSE": rmse, "ACC": acc})

        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(self.output_dir, outfile), index=False)
        return df

    @staticmethod
    def _mirror_merge(change_df: pd.DataFrame, pred_df: pd.DataFrame) -> pd.DataFrame:
        merged = pred_df.merge(change_df, left_on='region', right_on='region', how='inner').copy()
        if merged.empty:
            return pd.DataFrame(columns=['region','Clone1','Clone2','Clone1_CNA','Clone2_CNA','Clone1_predict_CNA','Clone2_predict_CNA'])
        def mode_for_row(row, prefix):
            cols = [c for c in merged.columns if c.startswith(row[prefix])]
            if not cols:
                return None
            vals = merged.loc[row.name, cols]
            m = vals.mode()
            return m.iloc[0] if not m.empty else None
        
        merged['Clone1_predict_CNA'] = merged.apply(mode_for_row, axis=1, prefix='Clone1')
        merged['Clone2_predict_CNA'] = merged.apply(mode_for_row, axis=1, prefix='Clone2')
        keep = ['region','Clone1','Clone2','Clone1_CNA','Clone2_CNA','Clone1_predict_CNA','Clone2_predict_CNA']
        return merged[keep].copy()

    def clusterConsistency(
        self,
        tool_clone_files: List[str],
        tool_names: List[str],
    ):
        result_list = []

        for path1, tool1 in zip(tool_clone_files, tool_names):

            if not os.path.exists(path1):
                result_list.append({
                    "Tool": tool1,
                    "ARI": None,
                    "AMI": None
                })
                continue

            data = pd.read_csv(path1)

            clusters1 = data['cell_id'].str.split("_").str[0]

            clusters2 = data['clone_id']

            # Compute the Adjusted Rand Index (ARI) between clusters and cell clones
            ari = adjustedRandIndex(clusters1, clusters2)
            
            # Compute the Adjusted Mutual Information (AMI) between clusters and cell clones
            ami = AMI(clusters1, clusters2)

                
            result_list.append({
                "Tool":tool1,
                "ARI": ari,
                "AMI": ami
            })

        result_df = pd.DataFrame(result_list)
        out = os.path.join(self.output_dir, "clustering_result.csv")
        result_df.to_csv(out,index=False)
        print(f"Clustering ARI,AMI saved to {out}") 

        return result_df

    def NA_ratio(
            self,
            tool_cna_files: List[str],
            tool_names: List[str],
            profile_bin_size = 100000,
            outfile: str = "NA_ratio_results.csv",
        ) -> pd.DataFrame:

        results = []

        for path, name in zip(tool_cna_files, tool_names):
            print(name)
            pred = pd.read_csv(path)

            pred = split_all_regions(pred.set_index("region"), profile_bin_size)
            pred = pred.reset_index().rename(columns={"index": "region"})

            na_ratio, na_cnt, total = evaluate_NA_ratio(pred)
            results.append({"Tool": name, "NA_ratio": na_ratio, "NA_count": na_cnt, "Total": total})
            print(results)

        df = pd.DataFrame(results)
        out = os.path.join(self.output_dir, outfile)
        df.to_csv(out, index=False)

        return df

    def detect_size(self,
        tool_cna_files: List[str],
        tool_names: List[str],
        size_file: str,
        level : str,
        outfile: str = "detect_size_acc.csv",
        profile_bin_size = 100000
        ):


        size_truth = read_and_drop_empty(size_file)

        size_truth = split_all_regions(size_truth.set_index("Segment"), profile_bin_size)
        size_truth = size_truth.reset_index().rename(columns={"index": "Segment"})

        results = []

        for f_h, name in zip(tool_cna_files, tool_names):            
            pred = read_and_drop_empty(f_h)

            pred = split_all_regions(pred.set_index("region"), profile_bin_size)
            pred = pred.reset_index().rename(columns={"index": "region"})

            if level == "cell":
                comb= self._size_join_cell(size_truth, pred)
            elif level == "clone":
                comb= self._size_join_clone(size_truth, pred)


            for t in comb['size'].unique().tolist():
                gt = comb[comb['size'] == t]['value']
                pd_ = comb[comb['size'] == t]['predict_value']
                acc = (pd_ == gt).mean()
                results.append({"Tool": name, "size": t, "ACC": acc})

        df = pd.DataFrame(results)
        df.to_csv(os.path.join(self.output_dir, outfile), index=False)
        return df

    @staticmethod        
    def _size_join_clone(change_df: pd.DataFrame, pred_df: pd.DataFrame) -> pd.DataFrame:
        
        merged = pred_df.merge(change_df, left_on='region', right_on='Segment', how='inner')
        if merged.empty:
            return pd.DataFrame(columns=['cell','Segment','value','predict_value','size'])

        cell_vals = merged['cell'].astype(str).to_numpy()

        all_cols = merged.columns.tolist()

        def build_prefix_map(prefixes):
            uniq = np.unique(prefixes)
            return {
                p: np.array([i for i, c in enumerate(all_cols) if c.startswith(p)], dtype=int)
                for p in uniq
            }

        cell_colmap = build_prefix_map(cell_vals)
        print("merged:" , merged.head())

        arr = merged.to_numpy()      # shape = (n_rows, n_cols)
        n = arr.shape[0]

        cell_pred = np.empty(n, dtype=object)

        for i in range(n):
            p = cell_vals[i]
            cols = cell_colmap.get(p)
            if cols is None or len(cols) == 0:
                cell_pred[i] = None
            else:
                cell_pred[i] = fast_mode_any(arr[i, cols])


        out = merged[[
            'cell','Segment','value','size'
        ]].copy()


        out['predict_value'] = cell_pred

        return out


    @staticmethod        
    def _size_join_cell(change_df: pd.DataFrame, pred_df: pd.DataFrame) -> pd.DataFrame:
        
        out = change_df[["cell", "Segment", "value", "size"]].copy()

        pred_wide = pred_df.set_index("region")  
        pred_cols = pred_wide.columns

        pred_value = np.full(len(out), None, dtype=object)

        cell_series = out["cell"].astype(str)
        seg_series = out["Segment"]

        for cell, idx in cell_series.groupby(cell_series).groups.items():
            if cell not in pred_cols:
                continue

            pred_value[idx] = pred_wide[cell].reindex(seg_series.iloc[idx]).to_numpy(dtype=object)

        out["predict_value"] = pred_value
        return out
    
    def cloneSizebycellprofile(
        self,
        gt_cna_file: str,
        tool_cna_files: List[str],
        tool_names: List[str],
        outfile: str = "clone_size_by_cell_profile.csv",
    ) -> pd.DataFrame:

        gt_df = pd.read_csv(gt_cna_file, index_col = 0)
        gt_df.fillna("1|1",inplace=True)

        gt_clone = get_cell_profile_size(gt_df)

        for f, name in zip(tool_cna_files, tool_names):
            pred_df = pd.read_csv(f,index_col = 0)
            pred_df = pred_df.dropna(axis=1, how='all')
            pred_df.fillna("1|1",inplace=True)

            pred_clone = get_cell_profile_size(pred_df)

            aligned = pred_clone.reindex(gt_clone.index)

            gt_clone[f"{name}_pred_size"] = aligned['cluster_size']
            # gt_clone[f"{name}_pred_cell"] = aligned.index

        gt_clone.to_csv(os.path.join(self.output_dir, outfile),index=True)

        size_col = "cluster_size"

        mean_cols = gt_clone.columns.drop(size_col)
        mean_cols = [c for c in mean_cols if pd.api.types.is_numeric_dtype(gt_clone[c])]

        df_mean = (
            gt_clone
            .groupby(size_col, as_index=False)[mean_cols]
            .mean()
            .rename(columns={size_col: "cluster_size"})
        )

        df_mean.to_csv(os.path.join(self.output_dir, f"mean_{outfile}"), index=False)

    def cloneSizebycluster(
        self,
        gt_cluster_file: str,
        tool_cluster_files: List[str],
        tool_names: List[str],
        outfile: str = "clone_size_by_cluster.csv",
    ) -> pd.DataFrame:

        gt_df = pd.read_csv(gt_cluster_file, index_col = 0)

        gt_clone = get_cluster_size(gt_df)

        for f, name in zip(tool_cluster_files, tool_names):
            pred_df = pd.read_csv(f,index_col = 0)

            pred_clone = get_cluster_size(pred_df)

            aligned = pred_clone.reindex(gt_clone.index)

            gt_clone[f"{name}_pred_size"] = aligned['cluster_size']
            # gt_clone[f"{name}_pred_cell"] = aligned.index

        gt_clone.to_csv(os.path.join(self.output_dir, outfile),index=True)

        size_col = "cluster_size"

        mean_cols = gt_clone.columns.drop(size_col)
        mean_cols = [c for c in mean_cols if pd.api.types.is_numeric_dtype(gt_clone[c])]

        df_mean = (
            gt_clone
            .groupby(size_col, as_index=False)[mean_cols]
            .mean()
            .rename(columns={size_col: "cluster_size"})
        )

        df_mean.to_csv(os.path.join(self.output_dir, f"mean_{outfile}"), index=False)

    # def segmentation(self,
    #     gt_cna_file: str,
    #     tool_cna_files: List[str],
    #     tool_names: List[str],
    #     profile_bin_size = 100000,
    #     threshold = 0.95,
    #     outprefix: str = "segmentation_",
    #     ):

    #     gt_profile =  pd.read_csv(gt_cna_file,index_col=0)
    #     gt_annotated_df = annotate_segments(gt_profile.T, detect_cnv_type=True, threshold=threshold)
    #     gt_annotated_df['region'] = gt_annotated_df['Chrom'] + ":" + gt_annotated_df['Start'].astype(str) + "-" + gt_annotated_df['End'].astype(str)

    #     gt_cna_event_counts = get_seg_cna_event_num(gt_annotated_df)
    #     gt_cna_event_counts = gt_cna_event_counts.rename(columns={"bin_level_count": f"gt_bin_level_count",
    #                                                               "seg_level_count": f"gt_seg_level_count"})


    #     overlap_results = []
    #     metric_results = []
    #     merged_cna_event_count = gt_cna_event_counts

        
    #     for f, name in zip(tool_cna_files, tool_names):

    #         tool_cna_df = pd.read_csv(f, index_col=0)
    #         tool_annotated_df = annotate_segments(tool_cna_df.T, detect_cnv_type=True, threshold=threshold)
    #         tool_annotated_df['region'] = tool_annotated_df['Chrom'] + ":" + tool_annotated_df['Start'].astype(str) + "-" + tool_annotated_df['End'].astype(str)
    #         print(tool_annotated_df.head())
    #         tool_annotated_df = split_all_regions(tool_annotated_df.set_index("region"), profile_bin_size)
    #         tool_annotated_df = tool_annotated_df.reset_index().rename(columns={"index": "region"})

    #         tool_counts = get_seg_cna_event_num(tool_annotated_df)
    #         tool_counts = tool_counts.rename(columns={"bin_level_count": f"{name}_bin_level_count",
    #                                                   "seg_level_count": f"{name}_seg_level_count"})

    #         merged_cna_event_count = merged_cna_event_count.merge(tool_counts, on=["size", "type"], how="outer")

    #         result = get_segment_overlap_ratio(gt_annotated_df, tool_annotated_df)
    #         result2 = get_segment_metric(gt_annotated_df,tool_annotated_df)
            

    #         result['Tool'] = name
    #         result2['Tool'] = name
    #         print(result)

    #         overlap_results.append(result)
    #         metric_results.append(result2)

    #     num_cols = [c for c in merged_cna_event_count.columns if c.endswith("_count")]
    #     merged_cna_event_count[num_cols] = merged_cna_event_count[num_cols].fillna(0).astype(int)

    #     overlap_df = pd.concat(overlap_results, ignore_index=True)
    #     metric_df = pd.concat(metric_results, ignore_index=True)
    #     overlap_df.to_csv(os.path.join(self.output_dir, f"{outprefix}overlap.csv"), index=False)
    #     metric_df.to_csv(os.path.join(self.output_dir, f"{outprefix}metrics.csv"), index=False)

    #     merged_cna_event_count.to_csv(os.path.join(self.output_dir, f"{outprefix}complex_cna_count.csv"), index=False)


    def segmentation(self,
        gt_cna_file: str,
        tool_cna_files: List[str],
        tool_names: List[str],
        profile_bin_size = 100000,
        threshold = 0.95,
        outprefix: str = "segmentation_",
        level = "cell",
        changes_file: Optional[str] = ''
        ):

        gt_profile =  pd.read_csv(gt_cna_file,index_col=0)
        gt_annotated_df = annotate_segments(gt_profile.T, detect_cnv_type=True, threshold=threshold)
        gt_annotated_df['region'] = gt_annotated_df['Chrom'] + ":" + gt_annotated_df['Start'].astype(str) + "-" + gt_annotated_df['End'].astype(str)

        metric_results = []
 
        for f, name in zip(tool_cna_files, tool_names):

            tool_cna_df = read_and_drop_empty(f)

            tool_cna_df_long = (
                pd.DataFrame(tool_cna_df)
                .melt(id_vars="region",
                        var_name="cell",
                        value_name="value")
                .dropna(subset=["value"])
                )

            tool_cna_df_long = split_all_regions(tool_cna_df_long.set_index("region"), profile_bin_size)
            tool_cna_df_long = tool_cna_df_long.reset_index().rename(columns={"index": "region"})
            print(tool_cna_df_long.head())

            result2 = get_segment_metric(gt_annotated_df,tool_cna_df_long)
            
            result2['Tool'] = name

            metric_results.append(result2)

        # if level == "clone":
        #     mirror_result = self.mirrorsubclone(
        #         tool_cna_files= tool_cna_files,
        #         tool_names= tool_names,
        #         changes_file= changes_file,
        #         profile_bin_size=profile_bin_size
        #     )
        #     mirror_result['Type'] = "Mirrored-Subclonal-CNA"

        #     metric_results.append(mirror_result)

        metric_df = pd.concat(metric_results, ignore_index=True)

        metric_df.to_csv(os.path.join(self.output_dir, f"{outprefix}metrics.csv"), index=False)


    # def seg_mirrored_subclonal(self,
    #     gt_cna_file: str,
    #     tool_cna_files: List[str],
    #     tool_names: List[str],
    #     outfile: str = "seg_mirrored_subclonal_count.csv",
    #     ):

    #     gt_profile =  read_and_drop_empty(gt_cna_file)
    #     gt_mirrored = find_mirrored_clones(gt_profile.set_index("region"))
        

    #     gt_mirrored_counts = get_seg_mirror_subclonal(gt_mirrored)
    #     gt_mirrored_counts = gt_mirrored_counts.rename(columns={"count": f"gt_count"})

    #     merged_mirrored_count = gt_mirrored_counts

    #     for f, name in zip(tool_cna_files, tool_names):

    #         tool_cna_df = read_and_drop_empty(f)
    #         tool_mirrored = find_mirrored_clones(tool_cna_df.set_index("region"))
        
    #         tool_mirrored_counts = get_seg_mirror_subclonal(tool_mirrored)
    #         tool_mirrored_counts = tool_mirrored_counts.rename(columns={"count": f"{name}_count"})

    #         merged_mirrored_count = merged_mirrored_count.merge(tool_mirrored_counts, on=["size"], how="outer")


    #     num_cols = [c for c in merged_mirrored_count.columns if c.endswith("_count")]
    #     merged_mirrored_count[num_cols] = merged_mirrored_count[num_cols].fillna(0).astype(int)
    #     merged_mirrored_count.to_csv(os.path.join(self.output_dir, f"{outfile}"), index=False)

    


    def rddetect(self,
        tool_cna_files: List[str],
        bin_count_files: List[str],
        tool_names: List[str],
        outfile: str = "rd_cn_l1.csv"
        ):

        results = {}

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

        return df


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