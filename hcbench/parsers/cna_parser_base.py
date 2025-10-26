import numpy as np
from .base import BaseParser
from .utils import split_and_process_haplotype
import pandas as pd
import os

class CNAParser(BaseParser):

    chrom_col: str
    start_col: str
    end_col: str
    cell_col: str
    value_col: str
    add_chr_prefix: bool = False
    start_plus_one: bool = False
    split_haplotype: bool = True   

    def preprocess_value(self, value):
        return value

    def before_pivot(self) -> pd.DataFrame:
        df = pd.read_csv(self.input_path, sep='\t')
        return df

    def run(self):

        df = self.before_pivot()

        if self.start_plus_one:
            df[self.start_col] = df[self.start_col].astype(int) + 1

        prefix = "chr" if self.add_chr_prefix else ""
        df['region'] = (
            prefix + df[self.chrom_col].astype(str)
            + ":" + df[self.start_col].astype(str)
            + "-" + df[self.end_col].astype(str)
        )
        

        df[self.value_col] = df[self.value_col].apply(self.preprocess_value)
        wide_df = df.pivot(index="region", columns=self.cell_col, values=self.value_col)

        wide_df['chr'] = wide_df.index.to_series().str.extract(r"chr?([0-9XY]+):")[0]
        wide_df['start'] = wide_df.index.to_series().str.extract(r":([0-9]+)-")[0].astype(int)
        wide_df.sort_values(['chr', 'start'], inplace=True)
        wide_df.drop(columns=['chr', 'start'], inplace=True)



        self._check_output_path()
        wide_df.to_csv(f"{self.output_path}/haplotype_combined.csv")
        print(f"[hcbench] {self.__class__.__name__} parsed CNA file saved to {self.output_path}/haplotype_combined.csv")

        if self.split_haplotype:
            self._postprocess_haplotype(wide_df)

    def _postprocess_haplotype(self, wide_df: pd.DataFrame):

        print(f"[hcbench] Splitting haplotypes → {self.output_path}")

        result = split_and_process_haplotype(wide_df)

        hap1 = result["hap1"]
        hap2 = result["hap2"]
        minor = result["minor"]
        major = result["major"]
        combined = result["combined"]

        hap1.to_csv(os.path.join(self.output_path, "haplotype_1.csv"))
        hap2.to_csv(os.path.join(self.output_path, "haplotype_2.csv"))

        minor.to_csv(os.path.join(self.output_path, "minor.csv"))
        major.to_csv(os.path.join(self.output_path, "major.csv"))
        combined.to_csv(os.path.join(self.output_path, "minor_major.csv"))

        print(f"[hcbench] ✅ Haplotype split complete. Files saved in {self.output_path}")
