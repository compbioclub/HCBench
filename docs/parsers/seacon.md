# SEACON Parser

This module provides `SeaconParser`, specialized for parsing SEACON outputs and exporting standardized matrices and helper files used by hcbench workflows.

Key features:

- Read SEACON CNA tables and automatically preprocess copy number values by replacing commas with pipes (`|`).
- generate `minor.csv`, `major.csv`, and combined `minor_major.csv` outputs by default.
- Export bin-level count matrices as `bin_counts.csv`, dynamically mapped to the standard regions extracted from the CNA output.
- Convert a headerless, tab-separated VAF long table into sparse matrix outputs.

## 🚀 Quick Start

### 🧬 1. Parse the CNA Matrix

```
from hcbench.parsers.seacon import SeaconParser

seacon_input = "/demo_output/seacon/cna_output.tsv"
seacon_output = "/output/seacon/"

seacon_parser = SeaconParser(
    input_path=seacon_input, 
    output_path=seacon_output,
)
seacon_parser.run()
```

After running, the parser will read the input file, split the haplotypes, and save the standardized files to the output directory. The output directory typically contains:

```
/output/seacon/
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
/output/seacon/
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
/output/seacon/
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
    chrom_col: str = "chrom",
    start_col: str = "start",
    end_col: str = "end",
    cell_col: str = "cell",
    value_col: str = "CN",
    start_offset: int = 0,
    add_chr_prefix: bool = False,
    output_haplotype: bool = False,
    **kwargs
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

























































## 📂 Input Files

The output directory of **SEACON** typically includes a single main file containing copy number information for all cells:

```
demo_output/seacon/
└── calls.tsv
```

- `calls.tsv` — the primary SEACON output file containing inferred copy number (CN) states per cell and genomic region.

An example of `calls.tsv`:

```
cell	chrom	start	end	CN
clone9_cell6	chr1	5000001	10000000	0,2
clone9_cell6	chr1	15000001	20000000	0,2
clone9_cell6	chr1	20000001	25000000	0,2
clone9_cell6	chr1	25000001	30000000	0,2
```

The **required columns** are:

```
chrom, start, end, cell, CN
```

Each entry in the `CN` column may contain **comma-separated allele-specific copy numbers** (e.g. `"0,2"`).
 The parser automatically converts these values to a standardized `"hap1|hap2"` format (e.g. `"0|2"`).

------

## 📤 Output Files

After parsing, the following file is generated in the specified output directory:

```
haplotype_combined.csv
```

- `haplotype_combined.csv` — standardized CNA matrix (regions × cells).
   Each value represents the haplotype-level copy number in the format `"hap1|hap2"`.

If `split_haplotype=True` is enabled in the base class,
 the parser will also produce additional derived matrices:

```
haplotype_1.csv
haplotype_2.csv
minor.csv
major.csv
minor_major.csv
```

------

## ⚙️ Key Parameters

| Parameter        | Description                                            | Default   |
| ---------------- | ------------------------------------------------------ | --------- |
| `chrom_col`      | Column name for chromosome.                            | `"chrom"` |
| `start_col`      | Column name for region start.                          | `"start"` |
| `end_col`        | Column name for region end.                            | `"end"`   |
| `cell_col`       | Column name for cell ID.                               | `"cell"`  |
| `value_col`      | Column name for copy number value.                     | `"CN"`    |
| `start_plus_one` | Whether to shift start coordinates by +1.              | `False`   |
| `add_chr_prefix` | Whether to enforce `"chr"` prefix on chromosome names. | `False`   |

------

## 🚀 Example Usage

```
from hcbench.parsers.seacon import SeaconParser

seacon_input = "/demo_output/seacon/calls.tsv"
seacon_output = "/output/seacon/"

seacon_parser = SeaconParser(
    input_path=seacon_input,
    output_path=seacon_output
)
seacon_parser.run()
```

After running, the parser will read `calls.tsv`,
 convert `CN` values from `"x,y"` to `"x|y"`,
 and save the standardized file:

```
haplotype_combined.csv
```

in your output directory.

------
