# hcbench/parsers/utils.py

import gzip
from pathlib import Path
import pandas as pd
import numpy as np
import os

def split_and_process_haplotype(wide_df: pd.DataFrame,change_hap : bool = False):

    hap1 = wide_df.applymap(lambda x: str(x).split('|')[0] if pd.notna(x) else None)
    hap2 = wide_df.applymap(lambda x: str(x).split('|')[1] if ('|' in str(x)) else None)

    hap1_num = hap1.apply(pd.to_numeric, errors="coerce")
    hap2_num = hap2.apply(pd.to_numeric, errors="coerce")

    hap1_num = hap1_num.fillna(-1)
    hap2_num = hap2_num.fillna(-1)

    minor = hap1_num.where(hap1_num <= hap2_num, hap2_num)
    major = hap1_num.where(hap1_num > hap2_num, hap2_num)

    minor = minor.replace(-1, np.nan).astype("Int64")
    major = major.replace(-1, np.nan).astype("Int64")

    combined = (
        minor.astype(str).replace("<NA>", "NA") + "|" +
        major.astype(str).replace("<NA>", "NA")
    )

    combined = combined.replace("NA|NA", "")


    return {
        "hap1": hap1,         
        "hap2": hap2,         
        "minor": minor,       
        "major": major,       
        "combined": combined  
    }

def split_region(row_label, bin_size):
        try:
            chrom, pos = row_label.split(":")
            start, end = map(int, pos.split("-"))
        except Exception:
            return [row_label]

        bins = []
        for s in range(start, end, bin_size):
            e = min(s + bin_size -1, end)
            bins.append(f"{chrom}:{s}-{e}")
        return bins

def split_all_regions(out_t, bin_size):

    labels = out_t.index.to_numpy()

    lens = np.fromiter(
        (len(split_region(l, bin_size)) for l in labels),
        dtype=int
    )
    
    new_index = []
    for l in labels:
        new_index.extend(split_region(l, bin_size))

    data = np.repeat(out_t.to_numpy(), lens, axis=0)

    out_new = pd.DataFrame(data, columns=out_t.columns, index=new_index)

    return out_new

def map_cell_to_barcode(df: pd.DataFrame, barcode_path: str, cell_col: str) -> pd.DataFrame:
    barcodes = pd.read_csv(barcode_path, sep='\t')


    required_cols = {"#CELL", "BARCODE"}
    if not required_cols.issubset(barcodes.columns):
        raise ValueError(f"❌ barcode.tsv must have {required_cols}")

    barcodes["#CELL"] = barcodes["#CELL"].astype(str).str.replace(r"\.bam$", "", regex=True)

    mapping = dict(zip(barcodes["BARCODE"].astype(str), barcodes["#CELL"].astype(str)))

    df[cell_col] = df[cell_col].astype(str).map(mapping).fillna(df[cell_col])

    matched = df[cell_col].isin(mapping.values()).sum()
    print(f"[hcbench] Mapped {matched} cells using barcode file.")
    
    return df


def read_table_auto(path: str) -> pd.DataFrame:
    """
    Read a CSV/TSV file with automatic delimiter detection.
    Supports .csv/.tsv/.txt and optionally gzip-compressed versions: .csv.gz/.tsv.gz/.txt.gz.
    """
    p = Path(path)
    suffixes = [s.lower() for s in p.suffixes]  # e.g. ['.csv', '.gz']

    is_gz = (len(suffixes) >= 1 and suffixes[-1] == ".gz")
    base_ext = suffixes[-2] if is_gz and len(suffixes) >= 2 else (suffixes[-1] if suffixes else "")

    # Choose separator based on base extension
    if base_ext == ".csv":
        sep = ","
    elif base_ext in [".tsv", ".txt"]:
        sep = "\t"
    else:
        # Auto-detect when extension is unknown
        open_fn = (lambda fp: gzip.open(fp, "rt")) if is_gz else (lambda fp: open(fp, "r", encoding="utf-8"))
        with open_fn(path) as f:
            first_line = f.readline()
        sep = "," if first_line.count(",") > first_line.count("\t") else "\t"
        print(f"Unknown file extension '{''.join(suffixes) or base_ext}', auto-detected separator as '{sep}'")

    # Read with selected separator; compression infer handles .gz automatically
    df = pd.read_csv(path, sep=sep, compression="infer")
    return df