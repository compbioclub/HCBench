import os
import numpy as np
import pandas as pd
from .cna_parser_base import CNAParser
from .utils import read_table_auto, split_all_regions
from ..utils import vcf_chr_files_to_ad_dp_mtx

class CNReinParser(CNAParser):
    chrom_col = "Chromosome"
    start_col = "Start"
    end_col = "End"
    cell_col = "Cell barcode"
    value_col = "HAP_CN"
    start_offset: int = 0
    add_chr_prefix = True


    def before_pivot(self):
        df = read_table_auto(self.input_path)
        df[self.value_col] = df["Haplotype 1"].astype(str) + "|" + df["Haplotype 2"].astype(str)

        cna_df = self.base_process(df)

        if self.bin_size is not None:
            cna_df = split_all_regions(cna_df, self.bin_size)

            cna_df.index.name = "region"

        return cna_df

    def get_bin_rdr(self, counts_path,cell_name_path, cna_path = None, 
                       ):
        
        if cna_path is None:
            cna_path = os.path.join(self.output_path, "haplotype_combined.csv")

        bin_cols = pd.read_csv(cna_path)['region'].tolist()

        counts_npz = np.load(counts_path)
        counts = counts_npz["arr_0"]  # shape: (n_cells, n_bins)

        cells_npz = np.load(cell_name_path)
        cell_names = cells_npz["arr_0"]  # shape: (n_cells,)

        print("Counts shape:", counts.shape)
        print("Cell names shape:", cell_names.shape)

        if counts.shape[0] != len(cell_names):
            raise ValueError(f"counts shape:  {counts.shape}, cell_names shape: {len(cell_names)}")

        # bin_cols = [f"bin_{i}" for i in range(counts.shape[1])]
        df = pd.DataFrame(counts.T, index=pd.Index(bin_cols, name="region"), columns=cell_names)

        output_file = os.path.join(self.output_path, "bin_rdr.csv")
        df.to_csv(output_file)
        print(f"[hcbench] Bin Rdr file saved to {output_file}")

    def get_VAF_matrix(self, vaf_file_dir, output_path = None,min_dp=1, min_cells=1,prefix = "cellSNP"):

        vcf_paths = [f"{vaf_file_dir}/seperates_chr{i}.vcf.gz" for i in range(1, 23)]
        missing = [p for p in vcf_paths if not os.path.exists(p)]
        if missing:
            raise FileNotFoundError(
                "Missing VCF files:\n" + "\n".join(missing)
            )
        
        if output_path is not None:
            output_vaf = os.path.join(output_path, "VAF")
        else:
            output_vaf = os.path.join(self.output_path, "VAF")

        vcf_chr_files_to_ad_dp_mtx(
            vcf_paths=vcf_paths,
            out_dir=output_vaf,
            prefix=f"{prefix}",
            keep_only_biallelic=True,
            min_dp=min_dp, min_cells=min_cells
        )


