# SEACON Parser

This module provides `SeaconParser`, specialized for parsing SEACON outputs and exporting standardized matrices and helper files used by hcbench workflows.

Key features:

- Read SEACON CNA tables and automatically preprocess copy number values by replacing commas with pipes (`|`).
- generate `minor.csv`, `major.csv`, and combined `minor_major.csv` outputs by default.
- Export bin-level count matrices as `bin_counts.csv`, dynamically mapped to the standard regions extracted from the CNA output.
- Convert a headerless, tab-separated VAF long table into sparse matrix outputs.

`bin_size` (positive integer, bp) and `genome` (`"hg19"` or `"hg38"`) are required. Creating the parser prints both values. `expand_bins=False` leaves input CNA regions intact and prints that they will not be split. Set `expand_bins=True` to split each input region into pieces of at most `bin_size` bp before writing CNA matrices; the parser prints that splitting is enabled. The pieces keep the input CNA values. `bin_size` also controls event-analysis resolution when splitting is off.

## 🚀 Quick Start

### 🧬 1. Parse the CNA Matrix

```
from pathlib import Path
from hcbench.parsers.seacon import SeaconParser

Path("output/seacon").mkdir(parents=True, exist_ok=True)

seacon_input = "/demo_output/seacon/calls.tsv"
seacon_output = "output/seacon/"

seacon_parser = SeaconParser(
    input_path=seacon_input, 
    output_path=seacon_output,
    bin_size=100000, genome="hg38",
)
seacon_parser.run()
```

After running, the parser will read the input file, split the haplotypes, and save the standardized files to the output directory. The output directory typically contains:

```
output/seacon/
├── minor.csv             
├── major.csv             
└── minor_major.csv       
```

- `minor_major.csv` — combined CNA matrix (regions × cells).

  Each value represents the combined haplotype copy number.

------

### 🧬 2. Parse the bin counts matrix

**⚠️ Important:** This method requires the `minor_major.csv` file to already exist in your output directory to properly extract the `region` index. Ensure you run the main parser pipeline first.

```
counts_file = "/demo_output/seacon/counts.tsv"

seacon_parser.get_bin_counts(counts_path=counts_file)
```

After running this command, the following standardized file will be created:

```
output/seacon/
└── bin_counts.csv
```

------

### 🧬 3. Parse the VAF sparse matrices

This method expects a tab-separated VAF file with **no header**.

Python

```
seacon_parser.get_VAF_matrix(
    vaf_file_path="/demo_output/seacon/vaf.tsv",
    min_dp=3,
    min_cells=10,
    prefix="cellSNP"
)
```

After running this command, the following standardized directory and Matrix Market files will be created:

Plaintext

```
output/seacon/
├── VAF/
│   ├── cellSNP_AD.mtx
│   ├── cellSNP_DP.mtx
│   └── ...
```

## ⚙️ Initialization

Python

```
SeaconParser(
    input_path: str,
    output_path: str,
    *, bin_size: int, genome: str,
    expand_bins: bool = False,
    output_haplotype: bool = False,
)
```

**`input_path`**: Path to the SEACON CNA output table.

The **required columns** based on the default configuration are:

```
chrom, start, end, cell, CN
```

**`output_path`**: Base output directory where the standardized matrices will be saved.

## 🧠 Core Methods

### `SeaconParser.run()`

Executes the standard pipeline for the CNA matrix:

------

### `SeaconParser.get_bin_counts(counts_path)`

Creates a region-by-cell wide matrix of per-bin counts.

**Input**

- `counts_path`: Path to the counts table. The parser expects a table with cells in rows, so it transposes the dataframe automatically.

**Output** writes to:

- `{self.output_path}/bin_counts.csv`.

This is a wide matrix:

- rows: `region` (dynamically matched from `minor_major.csv`).
- columns: cells.
- values: counts.

------

### `SeaconParser.get_VAF_matrix(vaf_file_path, output_path=None, min_dp=1, min_cells=1, prefix="cellSNP")`

Converts a tab-separated VAF table into sparse matrix outputs using `long_to_mtx()`.

**Input format**

`vaf_file_path` must be a **tab-separated file with no header** and exactly 5 columns. The parser automatically maps them to:

| **column index** | **meaning** |
| ---------------- | ----------- |
| 0                | `chr`       |
| 1                | `position`  |
| 2                | `cell`      |
| 3                | `Acount`    |
| 4                | `Bcount`    |

**Parameters**

`output_path` (optional)

- If provided: outputs under `{output_path}/VAF`, else outputs under `{self.output_path}/VAF`.

```
min_dp
```

- Filter low depth sites.

```
min_cells
```

- Filter sites supported by too few cells.

```
prefix
```

- Output file prefix (default: `cellSNP`).

**Output**

Creates a `VAF/` directory containing the Matrix Market (`.mtx`) files.
