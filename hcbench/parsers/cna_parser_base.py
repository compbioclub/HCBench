import numpy as np
from .base import BaseParser
from .utils import read_table_auto, split_all_regions, split_and_process_haplotype
import pandas as pd
import os

class CNAParser(BaseParser):

    chrom_col: str
    start_col: str
    end_col: str
    cell_col: str
    value_col: str
    add_chr_prefix: bool = False
    start_offset: int = 0
    split_haplotype: bool = True   


    def __init__(self, *args, ref_genome: str = "hg38",cell_list = None, bin_size = None,expand_bins: bool = False, **kwargs):
        self.ref_genome = ref_genome
        self.cell_list = cell_list
        self.bin_size = bin_size
        self.expand_bins = expand_bins
        super().__init__(*args, **kwargs)

    def preprocess_value(self, value):
        return value

    def _add_region_column(self, df: pd.DataFrame) -> pd.DataFrame:

        df[self.start_col] = df[self.start_col].astype(int) + self.start_offset

        prefix = "chr" if self.add_chr_prefix else ""
        df["region"] = (
            prefix + df[self.chrom_col].astype(str)
            + ":" + df[self.start_col].astype(str)
            + "-" + df[self.end_col].astype(str)
        )
        return df
    
    def _pivot_long_to_wide(self, df: pd.DataFrame) -> pd.DataFrame:
        df[self.value_col] = df[self.value_col].apply(self.preprocess_value)
        wide_df = df.pivot(index="region", columns=self.cell_col, values=self.value_col)
        return wide_df

    def _sort_wide_by_region(self, wide_df: pd.DataFrame) -> pd.DataFrame:
        wide_df['chr'] = wide_df.index.to_series().str.extract(r"chr?([0-9]+):")[0].astype(int)
        wide_df['start'] = wide_df.index.to_series().str.extract(r":([0-9]+)-")[0].astype(int)
        wide_df.sort_values(['chr', 'start'], inplace=True)
        wide_df.drop(columns=['chr', 'start'], inplace=True)

        return wide_df
    
    
    def base_process(self, df: pd.DataFrame) -> pd.DataFrame:
        df = self._add_region_column(df)
        wide_df = self._pivot_long_to_wide(df)
        wide_df = self._sort_wide_by_region(wide_df)
        return wide_df

    def before_pivot(self) -> pd.DataFrame:
        df = read_table_auto(self.input_path)
        cna_df = self.base_process(df)

        return cna_df

    def run(self):
        df = self.before_pivot()

        if self.expand_bins:
            df = self._build_cna_df(df)

        if self.cell_list is not None:
            df = df.reindex(columns=self.cell_list)
            # df = df.fillna("NA|NA")


        self._check_output_path()
        
        if self.split_haplotype:
            self._postprocess_haplotype(df)

        # df= df.replace("NA|NA", "")
        # df.to_csv(f"{self.output_path}/haplotype_combined.csv")
        print(f"[hcbench] {self.__class__.__name__} parsed CNA file saved to {self.output_path}/haplotype_combined.csv")


    def _postprocess_haplotype(self, wide_df: pd.DataFrame):

        print(f"[hcbench] Splitting haplotypes → {self.output_path}")

        result = split_and_process_haplotype(wide_df)

        if self.change_hap:
            hap1 = result["hap2"]
            hap2 = result["hap1"]
            hap_combined = (
                hap1.astype(str).replace("<NA>", "NA") + "|" +
                hap2.astype(str).replace("<NA>", "NA")
            )
            hap_combined= hap_combined.replace({"None|None": "", "nan|nan": "", "<NA>|<NA>": ""})

            hap_combined.to_csv(os.path.join(self.output_path, "haplotype_combined.csv"))

        else:
            hap1 = result["hap1"]
            hap2 = result["hap2"]

            hap_combined = (
                hap1.astype(str).replace("<NA>", "NA") + "|" +
                hap2.astype(str).replace("<NA>", "NA")
            )
            hap_combined= hap_combined.replace({"None|None": "", "nan|nan": "", "<NA>|<NA>": ""})

            hap_combined.to_csv(os.path.join(self.output_path, "haplotype_combined.csv"))
            
        minor = result["minor"]
        major = result["major"]
        combined = result["combined"]

        hap1.to_csv(os.path.join(self.output_path, "haplotype_1.csv"))
        hap2.to_csv(os.path.join(self.output_path, "haplotype_2.csv"))

        minor.to_csv(os.path.join(self.output_path, "minor.csv"))
        major.to_csv(os.path.join(self.output_path, "major.csv"))
        combined.to_csv(os.path.join(self.output_path, "minor_major.csv"))

        print(f"[hcbench] ✅ Haplotype split complete. Files saved in {self.output_path}")

    @staticmethod
    def _get_ref_path(filename):
        return os.path.join(os.path.dirname(os.path.dirname(__file__)), "ref", filename)

    def _build_cna_df(self,cna_df):
        print(f"[hcbench] Loading {self.ref_genome} centromere data")
        df_ref = pd.read_csv(self._get_ref_path(f"{self.ref_genome}_centromeres.csv"))

        df_ref["chrom"] = df_ref["chrom"].astype(str)

        df_ref["chrom_num"] = df_ref["chrom"].str.extract(r"(\d+)")

        df_ref["chrom_num"] = df_ref["chrom_num"].astype("float")  
        df_ref = df_ref[df_ref["chrom_num"].between(1, 22)].copy()


        df_ref['start'] = 1
        df_ref.rename(columns={"chrom_length": "end"},inplace=True)
        df_ref = df_ref[['chrom', 'start', 'end']]
        df_ref["region"] = df_ref.apply(
            lambda r: f"{r.chrom}:{r.start}-{r.end}", axis=1
        )

        df_bins = split_all_regions(df_ref.set_index("region"), self.bin_size)

        df_bins.index.name = "region"

        all_regions = df_bins.index

        cna_df_aligned = cna_df.reindex(all_regions)

        # cna_df_aligned = cna_df_aligned.fillna("NA|NA")

        return cna_df_aligned
