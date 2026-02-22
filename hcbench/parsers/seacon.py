import os
from .utils import read_table_auto, split_and_process_haplotype
from .cna_parser_base import CNAParser
import pandas as pd
from ..utils import long_to_mtx

class SeaconParser(CNAParser):
    chrom_col = "chrom"
    start_col = "start"
    end_col = "end"
    cell_col = "cell"
    value_col = "CN"
    start_offset: int = 0
    add_chr_prefix = False

    def __init__(self, *args, output_haplotype: bool = False, **kwargs):
        super().__init__(*args, **kwargs)
        self.output_haplotype = output_haplotype

    def preprocess_value(self, value):
        if isinstance(value, str):
            return value.replace(",", "|")
        return value

    def get_bin_counts(self, counts_path):

        df = read_table_auto(counts_path)

        df.set_index("cell", inplace=True)

        df = df.T

        df_cna= read_table_auto(path = os.path.join(self.output_path, "minor_major.csv"))

        region = df_cna["region"]

        df.index = region.tolist()

        df.index.name = "region"

        self._check_output_path()

        output_file = os.path.join(self.output_path, "bin_counts.csv")
        df.to_csv(output_file)

        print(f"[hcbench] Bin Counts file saved to {output_file}")
    
    def _postprocess_haplotype(self, wide_df: pd.DataFrame):

        if self.output_haplotype:
            return super()._postprocess_haplotype(wide_df)

        print(f"[hcbench] Splitting haplotypes → {self.output_path}")

        result = split_and_process_haplotype(wide_df)
            
        minor = result["minor"]
        major = result["major"]
        combined = result["combined"]

        # hap1.to_csv(os.path.join(self.output_path, "haplotype_1.csv"))
        # hap2.to_csv(os.path.join(self.output_path, "haplotype_2.csv"))

        minor.to_csv(os.path.join(self.output_path, "minor.csv"))
        major.to_csv(os.path.join(self.output_path, "major.csv"))
        combined.to_csv(os.path.join(self.output_path, "minor_major.csv"))

        print(f"[hcbench] split complete. Files saved in {self.output_path}")

    def get_VAF_matrix(self, vaf_file_path, output_path = None,min_dp=1, min_cells=1,prefix = "cellSNP"):
        
        df = pd.read_csv(vaf_file_path, sep="\t",header=None)

        df = df.rename(columns={
            0: "chr",
            1: "position",
            2: "cell",
            3: "Acount",
            4: "Bcount",
        })

        if output_path is not None:
            output_vaf = os.path.join(output_path, "VAF")
        else:
            output_vaf = os.path.join(self.output_path, "VAF")
        long_to_mtx(df,out_dir=output_vaf,min_dp=min_dp, min_cells=min_cells, prefix=prefix)




