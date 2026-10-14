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



    def annotate_segments(
        self, cna_df_path, threshold=0.8, detect_cnv_type=False,
        out_file="ground_truth_classified.csv", *, small_size=3_000_000,
        mid_size=10_000_000, chromosomes=None,
    ):
        """Keep the GT call defaults/path while using the shared annotation core."""
        super().annotate_segments(
            cna_df_path, threshold, detect_cnv_type, out_file,
            small_size=small_size, mid_size=mid_size,
            chromosomes=chromosomes,
        )

    # def segmentation(self,
    #     gt_cna_file: str,
    #     threshold = 0.8,
    #     outprefix: str = "segmentation_",
    #     ):

    #     gt_profile =  pd.read_csv(gt_cna_file,index_col=0)
    #     gt_annotated_df = annotate_segments(gt_profile.T, detect_cnv_type=True, threshold=threshold, genome=self.genome, bin_size=self.bin_size)
    #     gt_annotated_df['region'] = gt_annotated_df['Chrom'] + ":" + gt_annotated_df['Start'].astype(str) + "-" + gt_annotated_df['End'].astype(str)

    #     gt_cna_event_counts = get_seg_cna_event_num(gt_annotated_df)
    #     gt_cna_event_counts = gt_cna_event_counts.rename(columns={"bin_level_count": f"gt_bin_level_count",
    #                                                               "seg_level_count": f"gt_seg_level_count"})

    #     merged_cna_event_count = gt_cna_event_counts
    #     merged_cna_event_count.to_csv(os.path.join(self.output_path, f"{outprefix}complex_cna_count.csv"), index=False)


