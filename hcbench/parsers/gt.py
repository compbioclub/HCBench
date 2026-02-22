import os
import pandas as pd
import numpy as np
from .cna_parser_base import CNAParser
from .utils import split_and_process_haplotype
from ..utils import annotate_segments, get_seg_cna_event_num
import pkgutil

class GTParser(CNAParser):

    chrom_col = "Chromosome"
    start_col = "Start"
    end_col = "End"
    cell_col = "cell"
    value_col = "value"
    add_chr_prefix = False
    start_offset: int = 0
    annotate_segment = False
    split_haplotype: bool = True  

    @staticmethod
    def _split_segment(segment_str):
        try:
            chrom, coords = segment_str.split(":")
            start, end = map(int, coords.split("-"))
            return chrom, start, end
        except Exception as e:
            raise ValueError(f"Invalid segment format: {segment_str}") from e

    @staticmethod
    def _combined_segment(chrom, start, end):
        return f"{chrom}:{start}-{end}"


    def preprocess_value(self, value):
        if isinstance(value, str):
            return value.replace(",", "|")
        return value
    

    def base_process(self, df: pd.DataFrame) -> pd.DataFrame:
        df = self._add_region_column(df)
        df.set_index("region",inplace=True)
        df.drop(columns=[self.chrom_col,self.start_col,self.end_col],inplace=True)
        wide_df = self._sort_wide_by_region(df)
        return wide_df


    # def run(self):
    #     print(f"[hcbench] Parsing ground truth CNV: {self.input_path}")

    #     wide_df = pd.read_csv(self.input_path, index_col=0)

    #     result = split_and_process_haplotype(wide_df)

    #     os.makedirs(self.output_path, exist_ok=True)
    #     result["hap1"].to_csv(os.path.join(self.output_path, "haplotype_1.csv"))
    #     result["hap2"].to_csv(os.path.join(self.output_path, "haplotype_2.csv"))
    #     result["minor"].to_csv(os.path.join(self.output_path, "minor.csv"))
    #     result["major"].to_csv(os.path.join(self.output_path, "major.csv"))
    #     result["combined"].to_csv(os.path.join(self.output_path, "minor_major.csv"))

    #     print(f"[hcbench] Saved haplotype split results → {self.output_path}")

    #     if self.annotate_segment:

    #         annotated_df = self._annotate_segments(wide_df)
    #         annotated_path = os.path.join(self.output_path, "ground_truth_classified.csv")
    #         annotated_df.to_csv(annotated_path, index=False)
    #         print(f"[hcbench] Saved size/type classification → {annotated_path}")

    #         mirrored_df = self.find_mirrored_clones(wide_df)
    #         mirrored_path = os.path.join(self.output_path, "mirrored_clones.csv")
    #         mirrored_df.to_csv(mirrored_path, index=False)
    #         print(f"[hcbench] Saved mirrored clone events → {mirrored_path}")

    #         return {
    #             "hap1": result["hap1"],
    #             "hap2": result["hap2"],
    #             "minor": result["minor"],
    #             "major": result["major"],
    #             "combined": result["combined"],
    #             "classified": annotated_df,
    #             "mirrored": mirrored_df
    #         }
        
    #     else:
    #         return {
    #             "hap1": result["hap1"],
    #             "hap2": result["hap2"],
    #             "minor": result["minor"],
    #             "major": result["major"],
    #             "combined": result["combined"],
    #         }

    def annotate_segments(self, cna_df_path,threshold = 0.9,detect_cnv_type=False,out_file= "ground_truth_classified.csv"):
        df = pd.read_csv(cna_df_path)
        df.set_index("region",inplace=True)
        annotated_df = annotate_segments(df.T, threshold=threshold, detect_cnv_type=detect_cnv_type)
        annotated_path = os.path.join(self.output_path, out_file)
        annotated_df.to_csv(annotated_path, index=False)
        print(f"[hcbench] Saved size/type classification → {annotated_path}")

    # def _annotate_segments(
    #     self,
    #     df: pd.DataFrame,
    #     small_size: int = 3_000_000,
    #     mid_size: int = 10_000_000,
    #     genome: str = "hg38",
    # ) -> pd.DataFrame:
    #     print("[hcbench] Annotating ground truth CNV by size/type (whole-arm first)...")

    #     centromeres_df = self._load_centromeres(genome)
    #     merged_list = []
    #     cur_bin = []
        

    #     for cell, row in df.iterrows():
    #         chrom = start = end = cur_val = None

    #         for key, val in row.items():
    #             chrom_new = key.split(":")[0]
    #             start_new = key.split("-")[0].split(":")[-1]
    #             end_new = key.split("-")[-1]

    #             if val == "1|1" or val == "nan|nan" :
    #                 if chrom is not None:
    #                     cur_df = pd.DataFrame(cur_bin)
    #                     cur_df['size'] = self._classify_size(start, end, small_size, mid_size)
    #                     cur_df['type'] = self._classify_whole_arm_events(cur_val,chrom,start,end,centromeres_df)
    #                     merged_list.append(cur_df)
    #                     cur_bin = []
    #                 chrom = start = end = cur_val = None
    #                 continue

    #             if chrom is None or chrom_new != chrom or val != cur_val:
    #                 if len(cur_bin) >0:
    #                     cur_df = pd.DataFrame(cur_bin)
    #                     cur_df['size'] = self._classify_size(start, end,small_size, mid_size)
    #                     cur_df['type'] = self._classify_whole_arm_events(cur_val,chrom,start,end,centromeres_df)
    #                     merged_list.append(cur_df)
    #                 # print(f"----{len(cur_bin)} >0 ,chrom, start, end, cur_val = {chrom_new},{ start_new}, {end_new}, {val} ")
    #                 chrom, start, end, cur_val = chrom_new, start_new, end_new, val
    #                 cur_bin = [{"cell": cell, "Segment": f"{chrom}:{start}-{end}", "value": cur_val}]

    #             else:
    #                 end = end_new
    #                 cur_bin.append({"cell": cell, "Segment": f"{chrom}:{start_new}-{end_new}", "value": cur_val})


    #         # Finalize last segment for each cell
    #         if chrom is not None:
    #             cur_df = pd.DataFrame(cur_bin)
    #             cur_df['size'] = self._classify_size(start, end, small_size, mid_size)
    #             cur_df['type'] = self._classify_whole_arm_events(cur_val,chrom,start,end,centromeres_df)
    #             merged_list.append(cur_df)
    #             cur_bin =[]

    #     return pd.concat(merged_list, ignore_index=True)



        
    # @staticmethod
    # def _get_ref_path(filename):
    #     return os.path.join(os.path.dirname(os.path.dirname(__file__)), "ref", filename)

    # def _load_centromeres(self, genome="hg38"):
    #     print(f"[hcbench] Loading {genome} centromere data")
    #     df = pd.read_csv(self._get_ref_path(f"{genome}_centromeres.csv"))
    #     # print(f"[hcbench] Loading centromere data from {data_path}...")
    #     return df

    # @staticmethod
    # def _classify_size(start, end, small_size, mid_size):
    #     size = int(end) - int(start) + 1
    #     if size < small_size:
    #         return "focal"
    #     elif size <= mid_size:
    #         return "medium"
    #     return "broad"

    # @staticmethod
    # def _classify_cnv(val: str) -> str:
    #     try:
    #         left, right = map(int, str(val).split("|"))
    #     except Exception:
    #         return "UNKNOWN"

    #     s = left + right
    #     if s == 0:
    #         return "DEL"
    #     elif s == 1:
    #         return "CNL_LOH"
    #     elif s == 2 and (left == 0 or right == 0):
    #         return "CNN_LOH"
    #     elif s > 2 and (left == 0 or right == 0):
    #         return "CNG_LOH"
    #     elif s > 2 and (left != 0 and right != 0):
    #         return "DUP"
    #     return "UNKNOWN"

    # def _classify_whole_arm_events(
    #     self, val: str, chrom: str, start: int, end: int, centromeres_df: pd.DataFrame
    # ) -> str:
    #     try:
    #         left, right = map(int, str(val).split("|"))
    #     except Exception:
    #         return "UNKNOWN"

    #     ch = str(chrom)
    #     if not ch.startswith("chr"):
    #         ch = "chr" + ch

    #     row = centromeres_df.loc[centromeres_df["chrom"] == ch]
    #     if row.empty:
    #         return self._classify_cnv(val)
    #     row = row.iloc[0]

    #     p_end = int(row["p_arm_end"])
    #     q_start = int(row["q_arm_start"])
    #     chr_len = int(row["chrom_length"])

    #     covers_p = (start <= 1) and (end >= p_end)
    #     covers_q = (start <= q_start) and (end >= chr_len)

    #     if (left == 0 or right == 0) and ((start <= 1) and (end >= chr_len)):
    #         return "WCL"  # Whole-arm loss
    #     if (left > 1 or right > 1) and ((start <= 1) and (end >= chr_len)):
    #         return "WCD"  # Whole-arm gain
        
    #     if (left == 0 or right == 0) and (covers_p or covers_q):
    #         return "ARM-DEL"  # Whole-arm loss
    #     if (left > 1 or right > 1) and (covers_p or covers_q):
    #         return "ARM-DUP"  # Whole-arm gain

    #     return self._classify_cnv(val)

    def segmentation(self,
        gt_cna_file: str,
        threshold = 0.95,
        outprefix: str = "segmentation_",
        ):

        gt_profile =  pd.read_csv(gt_cna_file,index_col=0)
        gt_annotated_df = annotate_segments(gt_profile.T, detect_cnv_type=True, threshold=threshold)
        gt_annotated_df['region'] = gt_annotated_df['Chrom'] + ":" + gt_annotated_df['Start'].astype(str) + "-" + gt_annotated_df['End'].astype(str)

        gt_cna_event_counts = get_seg_cna_event_num(gt_annotated_df)
        gt_cna_event_counts = gt_cna_event_counts.rename(columns={"bin_level_count": f"gt_bin_level_count",
                                                                  "seg_level_count": f"gt_seg_level_count"})

        merged_cna_event_count = gt_cna_event_counts
        merged_cna_event_count.to_csv(os.path.join(self.output_path, f"{outprefix}complex_cna_count.csv"), index=False)


