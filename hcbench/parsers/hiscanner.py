import os

import pandas as pd
from .cna_parser_base import CNAParser
from functools import reduce

class HiscannerParser(CNAParser):
    """
    Parser for HiScanner output files.

    This class reads HiScanner-generated CNA files and converts them
    into a unified haplotype-level CNA table suitable for downstream analysis.
    """


    def __init__(
        self,
        hiscanner_final_cna_dir,
        **kwargs
    ):
        super().__init__(**kwargs)
        self.hiscanner_final_cna_dir = hiscanner_final_cna_dir

    def before_pivot(self):

        cell_list = [
            name.split('_input_table')[0]
            for name in os.listdir(self.hiscanner_final_cna_dir)
            if name.endswith("input_table.txt")
        ]

        dfs = []   

        for cell in cell_list:
            input_table_path = os.path.join(self.hiscanner_final_cna_dir, f"{cell}_input_table.txt")
            input_df = pd.read_csv(input_table_path, sep='\t')

            input_df['region'] = (
                input_df['CHROM'].astype(str)
                + ":" + input_df['START'].astype(str)
                + "-" + input_df['END'].astype(str)
            )

            input_df[cell] = input_df['A'].astype(str) + "|" + input_df['B'].astype(str)

            small_df = input_df[['region', cell]]

            dfs.append(small_df)

        merged_df = reduce(lambda left, right: pd.merge(left, right, on='region', how='outer'), dfs)
        merged_df.set_index('region', inplace=True)

        return merged_df


    def get_bin_counts(self, counts_col = "OBS"):

        cell_list = [
            name.split('_input_table')[0]
            for name in os.listdir(self.hiscanner_final_cna_dir)
            if name.endswith("input_table.txt")
        ]

        dfs = []   

        for cell in cell_list:
            input_table_path = os.path.join(self.hiscanner_final_cna_dir, f"{cell}_input_table.txt")
            input_df = pd.read_csv(input_table_path, sep='\t')

            input_df['region'] = (
                input_df['CHROM'].astype(str)
                + ":" + input_df['START'].astype(str)
                + "-" + input_df['END'].astype(str)
            )

            input_df[cell] = input_df[counts_col]

            small_df = input_df[['region', cell]]

            dfs.append(small_df)

        merged_df = reduce(lambda left, right: pd.merge(left, right, on='region', how='outer'), dfs)
        merged_df.set_index('region', inplace=True)

        self._check_output_path()
        output_file = os.path.join(self.output_path, "bin_counts.csv")
        merged_df.to_csv(output_file)
        print(f"[hcbench] Bin Counts file saved to {output_file}")
