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

    def annotate_segments(
        self, cna_df_path, threshold=0.9, detect_cnv_type=False,
        out_file="ground_truth_classified.csv", *, small_size=3_000_000,
        mid_size=10_000_000, genome=None, chromosomes=None,
    ):
        """Keep the GT call defaults/path while using the shared annotation core."""
        super().annotate_segments(
            cna_df_path, threshold, detect_cnv_type, out_file,
            small_size=small_size, mid_size=mid_size, genome=genome,
            chromosomes=chromosomes,
        )

    def segmentation(self,
        gt_cna_file: str,
        threshold = 0.95,
        outprefix: str = "segmentation_",
        ):

        gt_profile =  pd.read_csv(gt_cna_file,index_col=0)
        gt_annotated_df = annotate_segments(gt_profile.T, detect_cnv_type=True, threshold=threshold, genome=self.ref_genome)
        gt_annotated_df['region'] = gt_annotated_df['Chrom'] + ":" + gt_annotated_df['Start'].astype(str) + "-" + gt_annotated_df['End'].astype(str)

        gt_cna_event_counts = get_seg_cna_event_num(gt_annotated_df)
        gt_cna_event_counts = gt_cna_event_counts.rename(columns={"bin_level_count": f"gt_bin_level_count",
                                                                  "seg_level_count": f"gt_seg_level_count"})

        merged_cna_event_count = gt_cna_event_counts
        merged_cna_event_count.to_csv(os.path.join(self.output_path, f"{outprefix}complex_cna_count.csv"), index=False)


