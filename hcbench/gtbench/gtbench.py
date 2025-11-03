# hcbench/gtbench/gtbench.py
import os
import ast
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


class GTBench:

    def __init__(self, output_dir: str = "./output"):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)


    @staticmethod
    def _align(df_pred: pd.DataFrame, df_true: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        common_idx = df_pred.index.intersection(df_true.index)
        common_col = df_pred.columns.intersection(df_true.columns)
        return df_pred.loc[common_idx, common_col], df_true.loc[common_idx, common_col]

    @staticmethod
    def _flatten_pair(a: pd.DataFrame, b: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
        xa, xb = a.values.flatten().astype(float), b.values.flatten().astype(float)
        mask = ~np.isnan(xa) & ~np.isnan(xb)
        return xa[mask], xb[mask]

    @staticmethod
    def _safe_rmse(a: pd.DataFrame, b: pd.DataFrame) -> float:
        xa, xb = GTBench._flatten_pair(a, b)
        if xa.size == 0:
            return float("nan")
        return float(np.sqrt(mean_squared_error(xb, xa)))

    @staticmethod
    def _safe_spearman(a: pd.DataFrame, b: pd.DataFrame) -> float:
        xa, xb = GTBench._flatten_pair(a, b)
        if xa.size == 0:
            return float("nan")
        r, _ = spearmanr(xb, xa)
        return float(r)

    @staticmethod
    def _safe_acc(a: pd.DataFrame, b: pd.DataFrame) -> float:
        va, vb = a.values, b.values
        mask = ~pd.isna(va) & ~pd.isna(vb)
        n = mask.sum()
        if n == 0:
            return float("nan")
        return float((va[mask] == vb[mask]).mean())

    @staticmethod
    def _to_numeric(df: pd.DataFrame) -> pd.DataFrame:
        return df.apply(pd.to_numeric, errors="coerce")

    @staticmethod
    def _sum_from_combined(df: pd.DataFrame) -> pd.DataFrame:

        s0 = df.astype(str).apply(lambda x: x.str.split("|").str[0]).apply(pd.to_numeric, errors="coerce")
        s1 = df.astype(str).apply(lambda x: x.str.split("|").str[1]).apply(pd.to_numeric, errors="coerce")
        return s0 + s1

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
        outfile: str = "bin_level_results.csv",
        index_col: Optional[int] = 0,
    ) -> pd.DataFrame:

        truth = pd.read_csv(cna_profile_file, index_col=index_col)
        results = []

        for path, name in zip(tool_cna_files, tool_names):
            pred = pd.read_csv(path, index_col=index_col)
            rmse, scc, acc = self._evaluate_haplotype_predictions(pred, truth, haplotype)
            results.append({"Tool": name, "RMSE": rmse, "SCC": scc, "ACC": acc})

        df = pd.DataFrame(results)
        out = os.path.join(self.output_dir, outfile)
        df.to_csv(out, index=False)
        return df

    def _evaluate_haplotype_predictions(
        self, pred_df: pd.DataFrame, truth_df: pd.DataFrame, haplotype: str
    ) -> Tuple[float, float, float]:

        pred_df, truth_df = self._align(pred_df, truth_df)

        if haplotype == "combined":

            pred_sum = self._sum_from_combined(pred_df)
            truth_sum = self._sum_from_combined(truth_df)
            rmse = self._safe_rmse(pred_sum, truth_sum)
            scc  = self._safe_spearman(pred_sum, truth_sum)

            acc  = self._safe_acc_combined(pred_df, truth_df)
        else:

            pred_num = self._to_numeric(pred_df)
            truth_num = self._to_numeric(truth_df)
            rmse = self._safe_rmse(pred_num, truth_num)
            scc  = self._safe_spearman(pred_num, truth_num)
            acc  = self._safe_acc(pred_num, truth_num)

        return rmse, scc, acc
    
    def _safe_acc_combined(self, pred_df: pd.DataFrame, truth_df: pd.DataFrame) -> float:
        a_str = pred_df.astype(str)
        b_str = truth_df.astype(str)

        mask = (~pred_df.isna() & ~truth_df.isna()) & \
            (~a_str.isin(["nan", "nan|nan"]) & ~b_str.isin(["nan", "nan|nan"]))
        n = int(mask.values.sum())
        if n == 0:
            return float("nan")

        return float((a_str[mask].values == b_str[mask].values).mean())

    def cnclass(
        self,
        tool_hap1_cna_files: List[str],
        tool_hap2_cna_files: List[str],
        tool_names: List[str],
        profile_hap1_cna_file: str,
        profile_hap2_cna_file: str,
        type: str = "hcCNA",
        outfile: str = "cnclass_results.csv",
    ) -> pd.DataFrame:

        if type not in ["acCNA", "hcCNA"]:
            raise ValueError("type must be 'acCNA' or 'hcCNA'")
        hap_list = ["minor", "major"] if type == "acCNA" else ["hap1", "hap2"]

        gt_h1 = pd.read_csv(profile_hap1_cna_file)
        gt_h2 = pd.read_csv(profile_hap2_cna_file)

        conditions = [
            (">=2", "CN_Gain"), ("=1", "CN_Neutral"), ("=0", "CN_Loss"),
            ("=2", "CN_equal_2"), ("=3", "CN_equal_3"), ("=4", "CN_equal_4"),
            ("=5", "CN_equal_5"), ("=6", "CN_equal_6"), ("=7", "CN_equal_7"),
            ("=8", "CN_equal_8"), (">=9", "CN_over_9"),
        ]

        all_rows = []
        for f_h1, f_h2, tool in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
            p_h1 = pd.read_csv(f_h1)
            p_h2 = pd.read_csv(f_h2)
            for cond, folder in conditions:

                save_dir = os.path.join(self.output_dir, tool, folder)
                os.makedirs(save_dir, exist_ok=True)
                self._categorize_and_save(gt_h1, p_h1, save_dir, tool, cond, hap_list[0])
                self._categorize_and_save(gt_h2, p_h2, save_dir, tool, cond, hap_list[1])

                res = self._process_folder_for_metrics(save_dir, tool, hap_list)
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
    def _cat_fn_from_cond(condition: str):
        if condition.startswith(">="):
            thr = float(condition[2:])
            return lambda v: int(float(v) >= thr)
        if condition.startswith("="):
            thr = float(condition[1:])
            return lambda v: int(float(v) == thr)
        raise ValueError(f"no support condition: {condition}")

    def _categorize_and_save(self, df_truth: pd.DataFrame, df_pred: pd.DataFrame,
                             out_dir: str, tool: str, condition: str, haplo_type: str):
        cat = self._cat_fn_from_cond(condition)
        gt_c = df_truth.copy()
        pd_c = df_pred.copy()
        for col in df_truth.columns[1:]:
            gt_c[col] = df_truth[col].apply(cat)
            pd_c[col] = df_pred[col].apply(cat)
        gt_c.to_csv(os.path.join(out_dir, f"ground_truth_{haplo_type}.csv"), index=False)
        pd_c.to_csv(os.path.join(out_dir, f"{tool}_predict_{haplo_type}.csv"), index=False)

    def _process_folder_for_metrics(self, folder: str, tool: str, hap_list: List[str]) -> Dict[str, Dict[str, Dict]]:
        results = {}
        for h in hap_list:
            gpath = os.path.join(folder, f"ground_truth_{h}.csv")
            ppath = os.path.join(folder, f"{tool}_predict_{h}.csv")
            if not (os.path.exists(gpath) and os.path.exists(ppath)):
                continue
            gdf, pdf = pd.read_csv(gpath), pd.read_csv(ppath)
            gdf, pdf = self._align(gdf, pdf)
            results[h] = self._metrics_per_clone(gdf, pdf)
        return results

    @staticmethod
    def _align_for_csv(gdf: pd.DataFrame, pdf: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        cols = gdf.columns.intersection(pdf.columns)

        return gdf.loc[:, cols], pdf.loc[:, cols]

    @staticmethod
    def _metrics_per_clone(gdf: pd.DataFrame, pdf: pd.DataFrame) -> Dict[str, Dict]:
        out = {}

        clones = pdf.columns[1:].str.split("_").str[0].unique()
        for clone in clones:
            gsub = gdf.loc[:, gdf.columns.str.startswith(clone)]
            psub = pdf.loc[:, pdf.columns.str.startswith(clone)]
            out[clone] = GTBench._calc_bin_metrics(gsub, psub)
        return out

    @staticmethod
    def _calc_bin_metrics(gdf: pd.DataFrame, pdf: pd.DataFrame) -> Dict[str, Optional[float]]:
        y_true = gdf.values.flatten()
        y_pred = pdf.values.flatten()
        metrics = {}
        metrics["Accuracy"]  = accuracy_score(y_true, y_pred)
        metrics["Precision"] = precision_score(y_true, y_pred, zero_division=0)
        metrics["Recall"]    = recall_score(y_true, y_pred, zero_division=0)
        metrics["F1"]        = f1_score(y_true, y_pred, zero_division=0)
        metrics["Brier"]     = brier_score_loss(y_true, y_pred)

        if len(set(y_true)) > 1 and len(set(y_pred)) > 1:
            labels = sorted(set(y_true).union(set(y_pred)))
            precision, recall, _ = precision_recall_curve(y_true, y_pred)
            metrics["AUROC"] = roc_auc_score(y_true, y_pred)
            metrics["AUPRC"] = auc(recall, precision)
            tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=labels).ravel()
            metrics["Specificity"] = tn / (tn + fp) if (tn + fp) else None
            metrics["PPV"] = metrics["Precision"]
            metrics["NPV"] = tn / (tn + fn) if (tn + fn) else None
            metrics["Kappa"]     = cohen_kappa_score(y_true, y_pred)

        else:
            metrics["AUROC"] = None
            metrics["AUPRC"] = None
            metrics["Specificity"] = None
            metrics["PPV"] = None
            metrics["NPV"] = None
            metrics["Kappa"]  = None

        return metrics 

    # ========= 3) hccnchange =========
    def hccnchange(
        self,
        tool_hap1_cna_files: List[str],
        tool_hap2_cna_files: List[str],
        tool_names: List[str],
        changes_file: str,
        outfile: str = "evolution_onset_CN_Change.csv",
    ) -> pd.DataFrame:

        change_truth = pd.read_csv(changes_file)
    
        if 'Haplotype' in change_truth.columns:
            change_truth['Haplotype'] = self._normalize_hap_label(change_truth['Haplotype'])

        results = []

        for f_h1, f_h2, name in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
            p_h1 = pd.read_csv(f_h1)
            p_h2 = pd.read_csv(f_h2)

            h1 = self._hc_change_calc(change_truth, p_h1, 'hap1')
            h2 = self._hc_change_calc(change_truth, p_h2, 'hap2')
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
    def _hc_change_calc(change_df: pd.DataFrame, pred_df: pd.DataFrame, hap: str) -> pd.DataFrame:
        filt = change_df[GTBench._normalize_hap_label(change_df['Haplotype']) == hap]
        merged = pred_df.merge(filt, left_on='region', right_on='Segment', how='inner')
        if merged.empty:
            return pd.DataFrame(columns=['Parent','Child','Haplotype','Type','Segment','Change','Parent_predict_num','Child_predict_num','Change_predict'])

        def mode_for_row(row, prefix):
            cols = [c for c in merged.columns if c.startswith(row[prefix])]
            if not cols:
                return None
            vals = merged.loc[row.name, cols]
            m = vals.mode()
            return m.iloc[0] if not m.empty else None

        merged['Parent_predict_num'] = merged.apply(mode_for_row, axis=1, prefix='Parent')
        merged['Child_predict_num']  = merged.apply(mode_for_row, axis=1, prefix='Child')

        keep = ['Parent','Child','Haplotype','Type','Segment','Change','Parent_predict_num','Child_predict_num']
        out = merged[keep].copy()
        out['Change_predict'] = (
            out['Parent_predict_num'].fillna(0).astype(int).astype(str) + "->" +
            out['Child_predict_num'].fillna(0).astype(int).astype(str)
        )
        return out

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
    ) -> pd.DataFrame:

        tree = Phylo.read(tree_file, "newick")
        change_df = pd.read_csv(changes_file)
        if 'Haplotype' in change_df.columns:
            change_df['Haplotype'] = self._normalize_hap_label(change_df['Haplotype'])
        change_df = self._add_check_list(tree, change_df)

        results = []
        for f_h1, f_h2, name in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
            p_h1 = pd.read_csv(f_h1)
            p_h2 = pd.read_csv(f_h2)

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
        outfile: str = "evolution_onset_acc.csv",
    ) -> pd.DataFrame:

        change_truth = pd.read_csv(changes_file)
        if 'Haplotype' in change_truth.columns:
            change_truth['Haplotype'] = self._normalize_hap_label(change_truth['Haplotype'])
        results = []

        for f_h1, f_h2, name in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
            p_h1 = pd.read_csv(f_h1); p_h2 = pd.read_csv(f_h2)
            h1 = self._onset_join(change_truth, p_h1, 'hap1')
            h2 = self._onset_join(change_truth, p_h2, 'hap2')
            comb = pd.concat([h1, h2], ignore_index=True)

            for t in comb['Type'].unique().tolist():
                gt = comb[comb['Type'] == t]['Change']
                pd_ = comb[comb['Type'] == t]['Change_predict']
                acc = float((pd_ == gt).mean())
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
        def mode_for_row(row, prefix):
            cols = [c for c in merged.columns if c.startswith(row[prefix])]
            if not cols:
                return None
            vals = merged.loc[row.name, cols]
            m = vals.mode()
            return m.iloc[0] if not m.empty else None
        merged['Parent_predict_num'] = merged.apply(mode_for_row, axis=1, prefix='Parent')
        merged['Child_predict_num']  = merged.apply(mode_for_row, axis=1, prefix='Child')
        keep = ['Parent','Child','Haplotype','Type','Segment','Change','Parent_predict_num','Child_predict_num']
        out = merged[keep].copy()
        out['Change_predict'] = (
            out['Parent_predict_num'].fillna(0).astype(int).astype(str) + "->" +
            out['Child_predict_num'].fillna(0).astype(int).astype(str)
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
    ) -> pd.DataFrame:

        change_truth = pd.read_csv(changes_file)
        if 'Haplotype' in change_truth.columns:
            change_truth['Haplotype'] = self._normalize_hap_label(change_truth['Haplotype'])
        results = []

        for f_h1, f_h2, name in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
            p_h1 = pd.read_csv(f_h1); p_h2 = pd.read_csv(f_h2)
            h1 = self._onset_join(change_truth, p_h1, 'hap1')
            h2 = self._onset_join(change_truth, p_h2, 'hap2')
            comb = pd.concat([h1, h2], ignore_index=True)

            for t in comb['Type'].unique().tolist():
                gt = comb[comb['Type'] == t]['Change']
                pd_ = comb[comb['Type'] == t]['Change_predict']
                acc = self._parent_onset_acc(gt, pd_)
                pr, _ = self._parent_child_rmse(gt, pd_)
                results.append({"Tool": name, "Type": t, "RMSE": pr, "ACC": acc})

        df = pd.DataFrame(results)
        df.to_csv(os.path.join(self.output_dir, outfile), index=False)
        return df

    @staticmethod
    def _parent_onset_acc(gt: pd.Series, pd_: pd.Series) -> float:
        p  = gt.str.split('->').str[0].apply(pd.to_numeric, errors='coerce')
        pp = pd_.str.split('->').str[0].apply(pd.to_numeric, errors='coerce')
        mask = ~p.isna() & ~pp.isna()
        if not mask.any():
            return float("nan")
        return float((p[mask] == pp[mask]).mean())

    @staticmethod
    def _parent_child_rmse(gt: pd.Series, pd_: pd.Series) -> Tuple[float, float]:
        p  = gt.str.split('->').str[0].apply(pd.to_numeric, errors='coerce')
        pp = pd_.str.split('->').str[0].apply(pd.to_numeric, errors='coerce')
        c  = gt.str.split('->').str[1].apply(pd.to_numeric, errors='coerce')
        cp = pd_.str.split('->').str[1].apply(pd.to_numeric, errors='coerce')
        m1 = ~p.isna() & ~pp.isna()
        m2 = ~c.isna() & ~cp.isna()
        pr = float(np.sqrt(mean_squared_error(p[m1], pp[m1]))) if m1.any() else float("nan")
        cr = float(np.sqrt(mean_squared_error(c[m2], cp[m2]))) if m2.any() else float("nan")
        return pr, cr

    # ========= 7) hcPhasing =========
    def hcPhasing(
        self,
        tool_hap1_cna_files: List[str],
        tool_hap2_cna_files: List[str],
        tool_names: List[str],
        ground_truth_hap1_file: str,
        ground_truth_hap2_file: str,
        outfile: str = "hcPhasing.csv",
        index_col: Optional[int] = 0,
    ) -> pd.DataFrame:

        g1 = pd.read_csv(ground_truth_hap1_file, index_col=index_col)
        g2 = pd.read_csv(ground_truth_hap2_file, index_col=index_col)
        g1_bin, g2_bin = self._phase_to_binary(g1, g2)

        rows = []
        for f_h1, f_h2, name in zip(tool_hap1_cna_files, tool_hap2_cna_files, tool_names):
            t1 = pd.read_csv(f_h1, index_col=index_col)
            t2 = pd.read_csv(f_h2, index_col=index_col)
            h1_bin, h2_bin = self._phase_to_binary(t1, t2)
            part = self._eval_mismatch_switch(g1_bin, g2_bin, h1_bin, h2_bin, name)
            rows.extend(part)

        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(self.output_dir, outfile), index=False)
        return df

    @staticmethod
    def _phase_to_binary(hap1_df: pd.DataFrame, hap2_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
       
        h1 = (hap1_df > hap2_df).astype(int)
        h1 = h1.mask(hap1_df == hap2_df, -1).mask(hap1_df.isna() | hap2_df.isna())
        h2 = (hap2_df > hap1_df).astype(int)
        h2 = h2.mask(hap1_df == hap2_df, -1).mask(hap1_df.isna() | hap2_df.isna())
        return h1, h2

    @staticmethod
    def _eval_mismatch_switch(g1: pd.DataFrame, g2: pd.DataFrame,
                              h1: pd.DataFrame, h2: pd.DataFrame, tool_name: str) -> List[Dict]:
        g1 = g1.dropna(axis=1, how='all'); g2 = g2.dropna(axis=1, how='all')
   
        idx = g1.index.intersection(h1.index)
        col = g1.columns.intersection(h1.columns)
        g1, g2, h1, h2 = g1.loc[idx, col], g2.loc[idx, col], h1.loc[idx, col], h2.loc[idx, col]

        prefixes = col.str.extract(r"(^[^_]+)_")[0].unique()
        results = []

        for pref in prefixes:
            grp_cols = [c for c in col if c.startswith(f"{pref}_")]
            g1_g, g2_g = g1[grp_cols], g2[grp_cols]
            h1_g,  h2_g = h1[grp_cols],  h2[grp_cols]

            # mismatch
            mask_1 = g1_g != -1; mask_2 = g2_g != -1
            mm = (h1_g != g1_g)[mask_1].sum().sum()
            pm = (h2_g != g2_g)[mask_2].sum().sum()
            total_m = mask_1.sum().sum()
            # switch error
            nrow, ncol = g1_g.shape
            total_sw = 0; sw_1 = 0; sw_2 = 0
            for r in range(nrow):
                for c in range(ncol - 1):
                    if g1_g.iat[r, c] == -1 or g1_g.iat[r, c+1] == -1:
                        continue
                    total_sw += 1
                    if g1_g.iat[r, c] != h1_g.iat[r, c] or g1_g.iat[r, c+1] != h1_g.iat[r, c+1]:
                        sw_1 += 1
                    if g2_g.iat[r, c] != h2_g.iat[r, c] or g2_g.iat[r, c+1] != h2_g.iat[r, c+1]:
                        sw_2 += 1

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
        outfile: str = "mirror_subclone_result.csv",
    ) -> pd.DataFrame:

        change_df = pd.read_csv(changes_file)
        # change_df['region'] = (
        #     change_df['Chromosome'].astype(str) + ":" +
        #     change_df['Start'].astype(str) + "-" +
        #     change_df['End'].astype(str)
        # )

        rows = []
        for path, name in zip(tool_cna_files, tool_names):
            pred = pd.read_csv(path)
            combined = self._mirror_merge(change_df, pred)
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
        merged = pred_df.merge(change_df, left_on='region', right_on='region', how='inner')
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
