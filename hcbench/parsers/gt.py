import os
import pandas as pd
import numpy as np
from .cna_parser_base import CNAParser
from .utils import split_and_process_haplotype
import pkgutil

class GTParser(CNAParser):

    chrom_col = "chrom"
    start_col = "start"
    end_col = "end"
    cell_col = "cell"
    value_col = "value"
    add_chr_prefix = False
    start_plus_one = False

    @staticmethod
    def _split_segment(segment_str):
        try:
            chrom, coords = segment_str.split(":")
            start, end = map(int, coords.split("-"))
            return chrom, start, end
        except Exception as e:
            raise ValueError(f"Invalid segment format: {segment_str}") from e

    @staticmethod
    def _combined_segment(chrom, start, end):
        return f"{chrom}:{start}-{end}"


    def preprocess_value(self, value):
        if isinstance(value, str):
            return value.replace(",", "|")
        return value


    def run(self):
        print(f"[hcbench] Parsing ground truth CNV: {self.input_path}")

        wide_df = pd.read_csv(self.input_path, index_col=0)

        result = split_and_process_haplotype(wide_df)

        os.makedirs(self.output_path, exist_ok=True)
        result["hap1"].to_csv(os.path.join(self.output_path, "haplotype_1.csv"))
        result["hap2"].to_csv(os.path.join(self.output_path, "haplotype_2.csv"))
        result["minor"].to_csv(os.path.join(self.output_path, "minor.csv"))
        result["major"].to_csv(os.path.join(self.output_path, "major.csv"))
        result["combined"].to_csv(os.path.join(self.output_path, "minor_major.csv"))

        print(f"[hcbench] Saved haplotype split results → {self.output_path}")

        annotated_df = self._annotate_segments(wide_df)
        annotated_path = os.path.join(self.output_path, "ground_truth_classified.csv")
        annotated_df.to_csv(annotated_path, index=False)
        print(f"[hcbench] Saved size/type classification → {annotated_path}")

        mirrored_df = self.find_mirrored_clones(wide_df)
        mirrored_path = os.path.join(self.output_path, "mirrored_clones.csv")
        mirrored_df.to_csv(mirrored_path, index=False)
        print(f"[hcbench] Saved mirrored clone events → {mirrored_path}")

        return {
            "hap1": result["hap1"],
            "hap2": result["hap2"],
            "minor": result["minor"],
            "major": result["major"],
            "combined": result["combined"],
            "classified": annotated_df,
            "mirrored": mirrored_df
        }

    def _annotate_segments(
        self,
        wide_df: pd.DataFrame,
        small_size: int = 5_000_000,
        mid_size: int = 20_000_000,
        genome: str = "hg38",
    ) -> pd.DataFrame:
        print("[hcbench] Annotating ground truth CNV by size/type (whole-arm first)...")

        centromeres_df = self._load_centromeres(genome)
        records = []

        for region, row in wide_df.iterrows():
            chrom, coords = region.split(":")
            start, end = map(int, coords.split("-"))

            for cell, val in row.items():
                if pd.isna(val):
                    continue
                sval = str(val)
                if sval in ("1|1", "nan|nan"):
                    continue

                size_label = self._classify_size(start, end, small_size, mid_size)
                type_label = self._classify_whole_arm_events(
                    sval, chrom, start, end, centromeres_df
                )
                records.append(
                    {
                        "cell": cell,
                        "region": region,
                        "value": sval,
                        "size": size_label,
                        "type": type_label,
                    }
                )

        return pd.DataFrame(records)

        
    @staticmethod
    def _get_ref_path(filename):
        return os.path.join(os.path.dirname(os.path.dirname(__file__)), "ref", filename)

    def _load_centromeres(self, genome="hg38"):
        print(f"[hcbench] Loading {genome} centromere data")
        df = pd.read_csv(self._get_ref_path(f"{genome}_centromeres.csv"))
        # print(f"[hcbench] Loading centromere data from {data_path}...")
        return df

    @staticmethod
    def _classify_size(start, end, small_size, mid_size):
        size = int(end) - int(start) + 1
        if size <= small_size:
            return "small"
        elif size <= mid_size:
            return "middle"
        return "large"

    @staticmethod
    def _classify_cnv(val: str) -> str:
        try:
            left, right = map(int, str(val).split("|"))
        except Exception:
            return "UNKNOWN"

        s = left + right
        if s == 0:
            return "DEL"
        elif s == 1:
            return "CNL_LOH"
        elif s == 2 and (left == 0 or right == 0):
            return "CNN_LOH"
        elif s > 2 and (left == 0 or right == 0):
            return "CNG_LOH"
        elif s > 2 and (left != 0 and right != 0):
            return "DUP"
        return "UNKNOWN"

    def _classify_whole_arm_events(
        self, val: str, chrom: str, start: int, end: int, centromeres_df: pd.DataFrame
    ) -> str:
        try:
            left, right = map(int, str(val).split("|"))
        except Exception:
            return "UNKNOWN"

        ch = str(chrom)
        if not ch.startswith("chr"):
            ch = "chr" + ch

        row = centromeres_df.loc[centromeres_df["chrom"] == ch]
        if row.empty:
            return self._classify_cnv(val)
        row = row.iloc[0]

        p_end = int(row["p_arm_end"])
        q_start = int(row["q_arm_start"])
        chr_len = int(row["chrom_length"])

        covers_p = (start <= 1) and (end >= p_end)
        covers_q = (start <= q_start) and (end >= chr_len)

        if (left == 0 or right == 0) and (covers_p or covers_q):
            return "WCL"  # Whole-arm loss
        if (left > 1 or right > 1) and (covers_p or covers_q):
            return "WGD"  # Whole-arm gain

        return self._classify_cnv(val)

    @staticmethod
    def find_mirrored_clones(df: pd.DataFrame) -> pd.DataFrame:
        print("[hcbench] Detecting mirrored clone CNAs...")
        clone_cols = df.columns
        mirrored_rows = []

        for seg, row in df.iterrows():
            for i in range(len(clone_cols)):
                for j in range(i + 1, len(clone_cols)):
                    c1, c2 = row[clone_cols[i]], row[clone_cols[j]]
                    if pd.isna(c1) or pd.isna(c2):
                        continue
                    try:
                        hap1 = tuple(map(int, c1.split("|")))
                        hap2 = tuple(map(int, c2.split("|")))
                    except ValueError:
                        continue
                    if hap1 == hap2[::-1] and hap1[0] != hap1[1]:
                        mirrored_rows.append(
                            {
                                "region": seg,
                                "Clone1": clone_cols[i],
                                "Clone2": clone_cols[j],
                                "Clone1_CNA": c1,
                                "Clone2_CNA": c2,
                            }
                        )

        mirrored_df = pd.DataFrame(mirrored_rows)
        print(f"[hcbench] Found {len(mirrored_df)} mirrored clone events.")
        return mirrored_df
