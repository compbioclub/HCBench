# CHISEL Parser

## 📂 Input Files

The output directory of CHISEL typically looks like this:

```
demo_output/chisel/
├── calls/
│   └── calls.tsv
└── clones/
    └── mapping.tsv
```

- `calls.tsv` — the main file used to generate the CNA matrix.  
- `mapping.tsv` — an optional file mapping cells to their inferred clones.

An example of `calls.tsv`:

```
#CHR	START	END	CELL	NORM_COUNT	COUNT	RDR	A_COUNT	B_COUNT	BAF	CLUSTER	HAP_CN	CORRECTED_HAP_CN
chr1	0	1000000	AAACAGGTACAT	16269	1590	0.7594	76	67	0.4685	55	1|1	1|1
chr1	0	1000000	AAATTTGCCTTA	16269	3003	1.3324	74	195	0.7249	55	6|2	6|2
chr1	0	1000000	AACACATCCATC	16269	1587	0.9759	78	78	0.5	55	1|1	1|1
```

The **required columns** are:

```
#CHR, START, END, CELL, CORRECTED_HAP_CN
```

Please ensure these column names are spelled exactly as shown.

To maintain a consistent coordinate convention,  
the parser automatically **increments the `START` value by 1** — converting 0-based to 1-based coordinates.

---

## 📤 Output Files

After parsing, results are saved to the output directory, typically containing the following files:

```
haplotype_combined.csv
clusters.csv          # only if get_cluster() is used
haplotype_1.csv       # if split_haplotype=True
haplotype_2.csv       # if split_haplotype=True
minor.csv             # if split_haplotype=True
major.csv             # if split_haplotype=True
minor_major.csv       # if split_haplotype=True
```

- `haplotype_combined.csv` — main CNA matrix (regions × cells).
   Each value represents the combined haplotype copy number in the form `"hap1|hap2"`.
- `clusters.csv` — cell-to-clone mapping table produced by the `get_cluster()` function.

------

## ⚙️ Key Parameters

| Parameter        | Description                                        | Default |
| ---------------- | -------------------------------------------------- | ------- |
| `input_path`     | Path to the `calls.tsv` file.                      | —       |
| `output_path`    | Directory to save parsed results.                  | —       |
| `add_chr_prefix` | Whether to prepend `"chr"` to chromosome names.    | `False` |
| `start_plus_one` | Whether to convert `START` to 1-based coordinates. | `True`  |

------

## 🚀 Example Usage

### Parse the CNA Matrix

```
from hcbench.parsers.chisel import ChiselParser

chisel_input = "/demo_output/chisel/calls/calls.tsv"
chisel_output = "/output/chisel/"

chisel_parser = ChiselParser(chisel_input, chisel_output)
chisel_parser.run()
```

After running, the parser will read `calls.tsv` and automatically generate the standardized output files listed above.

------

### 🧬 Parse the Cluster File (New Feature)

If you have a CHISEL cluster mapping file (commonly named `mapping.tsv`),
 you can parse it separately using the new `get_cluster()` method:

```
cluster_file = "/demo_output/chisel/clones/mapping.tsv"
chisel_parser.get_cluster(cluster_file)
```

An example of `mapping.tsv`:

```
#CELL    CLUSTER    CLONE
AAACCTGAGAAGGACA    3233    Clone3233
AAACCTGAGATCTGCT    1484    Clone1484
AAACCTGAGTAATCCC    1924    None
AAACCTGAGTGCTGCC    3233    Clone3233

```

After running this command, the following standardized file will be created:

```
clusters.csv
```

with content:

```
cell_id,clone_id
AAACCTGAGAAGGACA,3233
AAACCTGAGATCTGCT,1484
AAACCTGAGTAATCCC,1924
```

This file ensures compatibility across tools in the HCBench framework.

------

## 🧠 Notes

- Missing columns in the cluster file (e.g., `#CELL` or `CLUSTER`) will trigger a clear error message.
- All outputs are written to the directory specified by `output_path`.
- If `split_haplotype=True`, additional per-haplotype and major/minor copy number matrices will also be generated.
