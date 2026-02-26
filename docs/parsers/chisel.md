# CHISEL Parser

This module provides ChiselParser, specialized for parsing CHISEL outputs and exporting standardized matrices and helper files used by hcbench workflows.

Key features:

- Run the existing  pipeline on CHISEL CNA tables
- Optional **two-pass** execution to export two different value columns into:
  - `clone_level/`
  - `cell_level/`
- Parse CHISEL cluster assignments into a standardized `clusters.csv`
- Export bin-level count matrices as `bin_counts.csv`
- Convert a cellSNP-like VAF long table into sparse matrix outputs

## 🚀 Quick Example

### 🧬 Parse the CNA Matrix

```
from hcbench.parsers.chisel import ChiselParser

chisel_input = "/demo_output/chisel/calls/calls.tsv"
chisel_output = "/output/chisel/"

chisel_parser = ChiselParser(chisel_input, chisel_output,value_cols=["HAP_CN_CLONE", "HAP_CN_CELL"])
chisel_parser.run()
```

After running, the parser will read `calls.tsv` and results are saved to the output directory, typically containing the following files:

```
/output/chisel/
├── clone-level/
│   └── haplotype_combined.csv
│   └── haplotype_1.csv       
│   └── haplotype_2.csv       
│   └── minor.csv             
│   └── major.csv             
│   └── minor_major.csv
├── cell-level/
│   └── haplotype_combined.csv
│   └── haplotype_1.csv       
│   └── haplotype_2.csv       
│   └── minor.csv             
│   └── major.csv             
│   └── minor_major.csv   
```

- `haplotype_combined.csv` — main CNA matrix (regions × cells).
  Each value represents the combined haplotype copy number in the form `"hap1|hap2"`.

------

### 🧬 Parse the Cluster File

If you have a CHISEL cluster mapping file (commonly named `mapping.tsv`), you can parse it separately using the new `get_cluster()` method:

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
/output/chisel/
└──clusters.csv
```

with content:

```
cell_id,clone_id
AAACCTGAGAAGGACA,3233
AAACCTGAGATCTGCT,1484
AAACCTGAGTAATCCC,1924
```

### 🧬 Parse the bin counts matrix

```
chisel_parser.get_bin_counts()
```

After running this command, the following standardized file will be created:

```
/output/chisel/
└──bin_counts.csv
```

### 🧬 Parse the VAF sparse matrices

```
chisel_parser.get_VAF_matrix(
    vaf_file_path="/demo_output/chisel/baf/baf.tsv",
    min_dp=3,
    min_cells=10,
)
```

After running this command, the following standardized file will be created:

```
/output/chisel/
├── VAF/
│   └──cellSNP_*.mtx
```



## ⚙️ Initialization

```python
ChiselParser(
    input_path: str,
    output_path: str,
    barcode_path: str | None = None,
    value_cols: list[str] | tuple[str, str] | None = None,
    **kwargs
)
```

**Parameters**

- **`input_path`**

   Path to the CHISEL CNA output table (long format), the output directory of CHISEL typically looks like this:

  ```
  demo_output/chisel/
  ├── calls/
  │   └── calls.tsv
  ```

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

​	Please ensure these column names are spelled exactly as shown.

- **`output_path`**
   Base output directory.
   If `value_cols` contains two columns, `ChiselParser` will automatically write to:

  - `{output_path}/clone_level`
  - `{output_path}/cell_level`

- **`barcode_path`** (optional)
   Path to a barcode mapping file used to remap cell.
   Mapping is applied in:

  - CNA table
  - cluster table
  - counts table

- **`value_cols`** (optional)
   Enables two-pass execution. Must be a `list`/`tuple` of **length 2**, e.g.:

```
  value_cols=['CORRECTED_HAP_CN','HAP_CN']
```

  Behavior:

  - `value_cols[0]` → written under `clone_level/`
  - `value_cols[1]` → written under `cell_level/`

  If not provided, the parser runs once using the default `value_col = "HAP_CN"`.



## 🧠 Core Methods

### `ChiselParser.run()`

- **Legacy mode (default)**: if `value_cols` is not provided (or only one column is configured), it runs once:
  - Uses `self.value_col` (default: `HAP_CN`)
  - Writes into `output_path`
- **Two-pass mode**: if `value_cols` has length 2, it runs twice:
  1. `self.value_col = value_cols[0]` → `output_path/clone_level`
  2. `self.value_col = value_cols[1]` → `output_path/cell_level`

------

### `ChiselParser.get_cluster(cluster_file_path)`

Parses a CHISEL cluster mapping file and writes a standardized CSV.

**Input**

- `cluster_file_path`: TSV file containing at least:
  - `#CELL`
  - `CLUSTER`

**Output**

Writes to:

- `{self.output_path}/clusters.csv`

CSV schema:

| column   | meaning                             |
| -------- | ----------------------------------- |
| cell_id  | cell identifier (possibly remapped) |
| clone_id | CHISEL cluster/clone label          |

------

### `ChiselParser.get_bin_counts()`

Creates a region-by-cell wide matrix of per-bin counts.

**Output**

Writes to:

- `{self.output_path}/bin_counts.csv`

This is a wide matrix:

- rows: `region` (e.g., `chr1:1000-2000`)
- columns: cells
- values: counts

------

### `ChiselParser.get_VAF_matrix(vaf_file_path, output_path=None, min_dp=1, min_cells=1, prefix="cellSNP")`

Converts a cellSNP-like VAF long table into sparse matrix outputs using `long_to_mtx()`.

**Input format**

`vaf_file_path` must be a **tab-separated file with no header** and exactly 5 columns:

| column index | meaning          |
| ------------ | ---------------- |
| 0            | chromosome       |
| 1            | genomic position |
| 2            | cell identifier  |
| 3            | allele A count   |
| 4            | allele B count   |

**Parameters**

- `output_path` (optional)
  - If provided: outputs under `{output_path}/VAF`
  - Else: outputs under `{self.output_path}/VAF`
- `min_dp`
  - filter low depth (typically `Acount + Bcount`)
- `min_cells`
  - filter sites supported by too few cells
- `prefix`
  - output file prefix (default: `cellSNP`)

**Output**

Creates a `VAF/` directory containing files:

  ```
.../VAF/
  cellSNP_*.mtx
  ...
  ```



