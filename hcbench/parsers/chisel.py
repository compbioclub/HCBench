import pandas as pd
from .cna_parser_base import CNAParser

class ChiselParser(CNAParser):
    chrom_col = "#CHR"
    start_col = "START"
    end_col = "END"
    cell_col = "CELL"
    value_col = "CORRECTED_HAP_CN"
    start_plus_one = True
    add_chr_prefix = False

    def get_cluster(self, cluster_file_path):
        print(f"[hcbench] Parsing Chisel cluster file...")
        
        cluster_df = pd.read_csv(cluster_file_path,sep='\t')

        result_df = pd.DataFrame({
            "cell_id": cluster_df['#CELL'],
            "clone_id": cluster_df['CLUSTER']
        })

        result_df.to_csv(f"{self.output_path}/clusters.csv", index=False)

    