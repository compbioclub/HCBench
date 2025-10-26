import pandas as pd
from .cna_parser_base import CNAParser

class CNReinParser(CNAParser):
    chrom_col = "Chromosome"
    start_col = "Start"
    end_col = "End"
    cell_col = "Cell barcode"
    value_col = "HAP_CN"
    start_plus_one = False
    add_chr_prefix = True


    def before_pivot(self):
        df = pd.read_csv(self.input_path)
        df[self.value_col] = df["Haplotype 1"].astype(str) + "|" + df["Haplotype 2"].astype(str)
        return df
