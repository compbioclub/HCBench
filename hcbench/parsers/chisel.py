import os
from pickle import FALSE
import pandas as pd
from .cna_parser_base import CNAParser
from .utils import map_cell_to_barcode, read_table_auto
from ..utils import long_to_mtx

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
    value_col = "HAP_CN"
    start_offset: int = 1
    add_chr_prefix = False

    def __init__(self, *args, barcode_path=None, value_cols=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.barcode_path = barcode_path

        # Backward-compatible behavior: if value_cols is not provided, run once with the default value_col
        if value_cols is None:
            self.value_cols = [self.value_col]
        else:
            if not isinstance(value_cols, (list, tuple)) or len(value_cols) != 2:
                raise ValueError("value_cols must be a list/tuple of length 2, e.g. ['colA', 'colB']")
            self.value_cols = list(value_cols)

        # Build two output directories from the given base output_path
        base_out = self.output_path
        self.output_paths = {
            "clone_level": f"{base_out}/clone_level",
            "cell_level": f"{base_out}/cell_level",
        }


    def before_pivot(self) -> pd.DataFrame:
        df = read_table_auto(self.input_path)

        if self.barcode_path is not None:
            df = map_cell_to_barcode(df, self.barcode_path, self.cell_col)

        cna_df = self.base_process(df)

        return cna_df

    def run(self):
        """
        Run the original CNAParser pipeline twice (without modifying CNAParser/BaseParser):
        - Once for the first value column -> {output_path}_clone_level
        - Once for the second value column -> {output_path}_cell_level

        Each run reuses CNAParser.run() exactly as-is; we only switch
        self.value_col and self.output_path before calling super().run().
        """
        # If only one column is configured (legacy mode), run the original pipeline once
        if len(self.value_cols) == 1:
            return super().run()

        # Fixed mapping per requirement: first column -> clone_level, second column -> cell_level
        mapping = [
            ("clone_level", self.value_cols[0]),
            ("cell_level", self.value_cols[1]),
        ]

        # Save the current state to avoid side effects if the same parser instance is reused
        old_value_col = self.value_col
        old_output_path = self.output_path

        try:
            for level_name, col in mapping:
                # Switch the value column and output path for this run
                self.value_col = col
                self.output_path = self.output_paths[level_name]
                if self.output_path and not os.path.exists(self.output_path):
                    os.makedirs(self.output_path, exist_ok=True)

                # Execute the original CNAParser logic: read -> pivot -> (optional) expand bins
                # -> (optional) reindex cells -> haplotype postprocess -> write outputs
                super().run()
        finally:
            # Restore state
            self.value_col = old_value_col
            self.output_path = old_output_path


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

        if self.barcode_path is not None:
            cluster_df = map_cell_to_barcode(cluster_df, self.barcode_path, "#CELL")


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

    def get_bin_counts(self, counts_col = "COUNT"):

        df = read_table_auto(self.input_path)

        if self.barcode_path is not None:
            df = map_cell_to_barcode(df, self.barcode_path, self.cell_col)

        df = self._add_region_column(df)

        wide_df = df.pivot(index="region", columns=self.cell_col, values=counts_col)

        wide_df['chr'] = wide_df.index.to_series().str.extract(r"chr?([0-9XY]+):")[0]
        wide_df['start'] = wide_df.index.to_series().str.extract(r":([0-9]+)-")[0].astype(int)
        wide_df.sort_values(['chr', 'start'], inplace=True)
        wide_df.drop(columns=['chr', 'start'], inplace=True)



        self._check_output_path()

        output_file = os.path.join(self.output_path, "bin_counts.csv")
        wide_df.to_csv(output_file)

        print(f"[hcbench] Bin Counts file saved to {output_file}")

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


