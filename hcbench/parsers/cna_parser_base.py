from .base import BaseParser
import pandas as pd

class CNAParser(BaseParser):

    chrom_col: str
    start_col: str
    end_col: str
    cell_col: str
    value_col: str
    add_chr_prefix: bool = False
    start_plus_one: bool = False

    def preprocess_value(self, value):
        return value

    def run(self):
        df = pd.read_csv(self.input_path, sep='\t')

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

        wide_df.to_csv(self.output_path)
        print(f"[hcbench] {self.__class__.__name__} parsed CNA file saved to {self.output_path}")
