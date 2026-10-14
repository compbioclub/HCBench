# SIGNALS Parser

This module provides `SignalsParser`, specialized for parsing Signals outputs and exporting standardized matrices and helper files used by hcbench workflows.

Key features:

- Run the existing pipeline on Signals CNA tables
- Parse Signals cluster assignments into a standardized `clusters.csv`
- Export bin-level count matrices as `bin_counts.csv`
- Convert a VAF long table into sparse matrix outputs

`bin_size` (positive integer, bp) and `genome` (`"hg19"` or `"hg38"`) are required. Creating the parser prints both values. `expand_bins=False` leaves input CNA regions intact and prints that they will not be split. Set `expand_bins=True` to split each input region into pieces of at most `bin_size` bp before writing CNA matrices; the parser prints that splitting is enabled. The pieces keep the input CNA values. `bin_size` also controls event-analysis resolution when splitting is off.

## 🚀 Quick Start

### 🧬 Exporting hscn_data.tsv from R

Run the following R code to extract and save the `hscn$data` component:

```
mnt <- "/demo_output/signals/output/"

hscn <- readRDS("/demo_output/signals/hscn.rds")

# Export the data component
write.table(
  hscn$data,
  file = paste0(mnt, "hscn_data.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
```

This will produce the file `hscn_data.tsv`, which is then used as input for the parser.

### 🧬 1. Parse the CNA Matrix

```python
from pathlib import Path
from hcbench.parsers.signals import SignalsParser

Path("output/signals").mkdir(parents=True, exist_ok=True)

signals_input = "/demo_output/signals/output/hscn_data.tsv"
signals_output = "output/signals/"

signals_parser = SignalsParser(input_path=signals_input, output_path=signals_output, bin_size=100000, genome="hg38")
signals_parser.run()
```

After running, the parser will read the input file and results are saved to the output directory, typically containing the following files:

```
output/signals/
├── haplotype_combined.csv
├── haplotype_1.csv       
├── haplotype_2.csv       
├── minor.csv             
├── major.csv             
└── minor_major.csv
```

- `haplotype_combined.csv` — main CNA matrix (regions × cells).



### 🧬 2. Parse the Cluster File

If you have a Signals cluster mapping file, you can parse it separately using the `get_cluster()` method:

```
cluster_file = "/demo_output/signals/clusters.csv"
signals_parser.get_cluster(cluster_file)
```

An example of the input cluster file:

```
cell_id,clone_id
AAACCTGAGAAGGACA,CloneA
AAACCTGAGATCTGCT,CloneB
AAACCTGAGTAATCCC,CloneA
```

After running this command, the following standardized file will be created:

```
output/signals/
└── clusters.csv
```



### 🧬 3. Parse the bin counts matrix

Unlike the CHISEL parser, the Signals parser requires the explicit path to the bin counts file:

```
bin_count_file = "hmmcopy_results/reads.csv.gz"
signals_parser.get_bin_counts(bin_count_file)
```

After running this command, the following standardized file will be created:

```
output/signals/
└── bin_counts.csv
```

------

### 🧬 4. Parse the VAF sparse matrices

```
signals_parser.get_VAF_matrix(
    vaf_file_path="hscn_pipeline_apptainer/results/counthaps/allele_counts_all.csv.gz",
    min_dp=3,
    min_cells=10,
)
```

After running this command, the following standardized directory and files will be created:

```
output/signals/
├── VAF/
│   └── cellSNP_*.mtx
```



## ⚙️ Initialization

```
SignalsParser(
    input_path: str,
    output_path: str,
    *, bin_size: int, genome: str,
    expand_bins: bool = False,
)
```

**`input_path`**: Path to the Signals CNA output. An example of the expected input format:

```
chr start   end reads   copy    state   cell_id alleleA alleleB totalcounts BAF state_min   A   B   state_AS_phased state_AS    LOH phase   state_phase state_BAF
1   5000001 10000000    34762   NA  2   clone1_cell1    677 193 870 0.22183908045977    1   1   1   1|1 1|1 NO  Balanced    Balanced    0.5
1   20000001    25000000    34200   NA  2   clone1_cell1    639 222 861 0.257839721254355   1   1   1   1|1 1|1 NO  Balanced    Balanced    0.5
1   30000001    35000000    42510   NA  2   clone1_cell1    807 203 1010    0.200990099009901   1   1   1   1|1 1|1 NO  Balanced    Balanced    0.5
```

The **required columns** are:

```
chr, start, end, cell_id, state_AS_phased
```

**`output_path`**: Base output directory where the standardized matrices will be saved.

## 🧠 Core Methods

### `SignalsParser.run()`

Executes the standard pipeline for the CNA matrix:

------

### `SignalsParser.get_cluster(cluster_file_path)`

Parses a Signals cluster mapping file and writes a standardized CSV.

**Input**

- `cluster_file_path`: CSV file containing at least the columns: `cell_id`, `clone_id`.

**Output** writes to:

- `{self.output_path}/clusters.csv`, strictly retaining these two columns.

------

### `SignalsParser.get_bin_counts(bin_count_file_path)`

Creates a region-by-cell wide matrix of per-bin counts. Note that this initializes an internal temporary parser specifically configured for the parser's `genome` and targets the `reads` column without adding chromosome prefixes.

**Input**

- `bin_count_file_path`: Path to the raw bin counts CSV file.

**Output** writes to:

- `{self.output_path}/bin_counts.csv`

This is a wide matrix:

- rows: `region`
- columns: cells
- values: counts

------

### `SignalsParser.get_VAF_matrix(vaf_file_path, output_path=None, min_dp=1, min_cells=1, prefix="cellSNP")`

Converts a VAF long table into sparse matrix outputs.

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

Creates a `VAF/` directory containing the Matrix Market (`.mtx`) files:

Plaintext

```
.../VAF/
└── cellSNP_*.mtx
```



## 🧩 Example of hscn_data.tsv

An example of the exported file may look like:

```
chr	start	end	reads	copy	state	cell_id	alleleA	alleleB	totalcounts	BAF	state_min	A	B	state_AS_phased	state_AS	LOH	phase	state_phase	state_BAF
1	5000001	10000000	34762	NA	2	clone1_cell1	677	193	870	0.22183908045977	1	1	1	1|1	1|1	NO	Balanced	Balanced	0.5
1	20000001	25000000	34200	NA	2	clone1_cell1	639	222	861	0.257839721254355	1	1	1	1|1	1|1	NO	Balanced	Balanced	0.5
1	30000001	35000000	42510	NA	2	clone1_cell1	807	203	1010	0.200990099009901	1	1	1	1|1	1|1	NO	Balanced	Balanced	0.5
```

Each row corresponds to one genomic bin in a given single cell.

The ***\*required columns\**** are:

\```

chr,start,end,cell_id,state_AS_phased

\```

