import os
import pandas as pd
from .cna_parser_base import CNAParser


class ChiselParser(CNAParser):
    """
    Parser for CHISEL output files.

    This class extends the generic CNAParser and provides a method
    to parse CHISEL cluster mapping files.
    """

    chrom_col = "#CHR"
    start_col = "START"
    end_col = "END"
    cell_col = "CELL"
    value_col = "CORRECTED_HAP_CN"
    start_plus_one = True
    add_chr_prefix = False

    def get_cluster(self, cluster_file_path):
        """
        Parse the CHISEL cluster mapping file and save a standardized CSV.

        Args:
            cluster_file_path (str): Path to the CHISEL cluster file.

        Output:
            clusters.csv — contains two columns:
                cell_id, clone_id
        """
        print("[hcbench] Parsing CHISEL cluster file...")

        try:
            cluster_df = pd.read_csv(cluster_file_path, sep="\t")
        except Exception as e:
            raise ValueError(f"Failed to read cluster file: {e}")

        required_cols = {"#CELL", "CLUSTER"}
        missing_cols = required_cols - set(cluster_df.columns)
        if missing_cols:
            raise ValueError(f"Missing required columns: {', '.join(missing_cols)}")

        result_df = pd.DataFrame({
            "cell_id": cluster_df["#CELL"],
            "clone_id": cluster_df["CLUSTER"]
        })

        output_file = os.path.join(self.output_path, "clusters.csv")
        result_df.to_csv(output_file, index=False)

        print(f"[hcbench] Cluster file saved to {output_file}")
