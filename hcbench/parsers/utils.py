# hcbench/parsers/utils.py

import pandas as pd
import numpy as np
import os

def split_and_process_haplotype(wide_df: pd.DataFrame):

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
            e = min(s + bin_size, end)
            bins.append(f"{chrom}:{s}-{e}")
        return bins