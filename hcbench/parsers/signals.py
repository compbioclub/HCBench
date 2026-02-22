import os
import pandas as pd
from .cna_parser_base import CNAParser
from ..utils import long_to_mtx

class SignalsParser(CNAParser):
    chrom_col = "chr"
    start_col = "start"
    end_col = "end"
    cell_col = "cell_id"
    value_col = "state_AS_phased"
    start_offset: int = 0
    add_chr_prefix = True


    def get_cluster(self, cluster_file_path):
        """
        Parse the Signal cluster mapping file and save a standardized CSV.

        Args:
            cluster_file_path (str): Path to the Signal cluster file.

        Output:
            clusters.csv — contains two columns:
                cell_id, clone_id
        """
        print("[hcbench] Parsing Signal cluster file...")

        try:
            cluster_df = pd.read_csv(cluster_file_path)
        except Exception as e:
            raise ValueError(f"Failed to read cluster file: {e}")

        required_cols = {"cell_id", "clone_id"}

        missing_cols = required_cols - set(cluster_df.columns)
        if missing_cols:
            raise ValueError(f"Missing required columns: {', '.join(missing_cols)}")

        result_df = pd.DataFrame({
            "cell_id": cluster_df["cell_id"],
            "clone_id": cluster_df["clone_id"]
        })

        output_file = os.path.join(self.output_path, "clusters.csv")
        result_df.to_csv(output_file, index=False)

        print(f"[hcbench] Cluster file saved to {output_file}")


    def get_bin_counts(self, bin_count_file_path):
        
        p = SignalsParser(ref_genome="hg38",cell_list=self.cell_list,input_path = bin_count_file_path, output_path = self.output_path)
        p.value_col = 'reads'
        p.add_chr_prefix = False
        df = p.before_pivot()

        output_file = os.path.join(self.output_path, "bin_counts.csv")
        df.to_csv(output_file)
        print(f"[hcbench] Bin Counts file saved to {output_file}")

    def get_VAF_matrix(self, vaf_file_path, output_path = None,min_dp=1, min_cells=1,prefix = "cellSNP"):
        
        df = pd.read_csv(vaf_file_path, index_col=0)

        if output_path is not None:
            output_vaf = os.path.join(output_path, "VAF")
        else:
            output_vaf = os.path.join(self.output_path, "VAF")
        long_to_mtx(df,out_dir=output_vaf,chr_col="chr",pos_col="position",cell_col="cell_id",a_col="allele0",b_col="allele1",ad_from="B",min_dp=min_dp, min_cells=min_cells, prefix=prefix)

