# SIGNALS Parser

This module provides `SignalsParser`, specialized for parsing Signals outputs and exporting standardized matrices and helper files used by hcbench workflows.



































## 📂 Input Files

The output directory of **SIGNALS** typically contains the R object `hscn.rds`, which stores haplotype-specific copy number data.

```
demo_output/signals/output/
├── hscn.rds
└── hscn_data.tsv   # exported for Python parser
```

- `hscn.rds` — main SIGNALS output file (R serialized object).
- `hscn_data.tsv` — exported tab-delimited file containing per-cell phased copy number states.

Because the SIGNALS parser in HCBench operates in Python,
 you need to **first export** the data component from the RDS file into a `.tsv` file using R.

------

## 🧬 Exporting hscn_data.tsv from R

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

------

## 📤 Output Files

After parsing, results are saved in the specified output directory:

```
haplotype_combined.csv
```

- `haplotype_combined.csv` — standardized CNA matrix (regions × cells).
   Each entry represents the phased haplotype-specific copy number state of a cell.

If `split_haplotype=True` is enabled in the base class,
 additional files such as `haplotype_1.csv`, `haplotype_2.csv`, and `minor_major.csv` will also be generated.

------

## ⚙️ Key Parameters

| Parameter        | Description                                            | Default             |
| ---------------- | ------------------------------------------------------ | ------------------- |
| `chrom_col`      | Column name for chromosome.                            | `"chr"`             |
| `start_col`      | Column name for region start.                          | `"start"`           |
| `end_col`        | Column name for region end.                            | `"end"`             |
| `cell_col`       | Column name for cell identifier.                       | `"cell_id"`         |
| `value_col`      | Column name for the haplotype-phased CN value.         | `"state_AS_phased"` |
| `start_plus_one` | Whether to convert 0-based to 1-based coordinates.     | `False`             |
| `add_chr_prefix` | Whether to ensure chromosome names begin with `"chr"`. | `True`              |

------

## 🚀 Example Usage

```
from hcbench.parsers.signals import SignalsParser

signals_input = "/demo_output/signals/output/hscn_data.tsv"
signals_output = "/output/signals/"

signals_parser = SignalsParser(
    input_path=signals_input,
    output_path=signals_output
)
signals_parser.run()
```

After running, the parser reads the exported `hscn_data.tsv`, extracts the `"state_AS_phased"` column as the haplotype-level copy number signal, and saves the standardized file:

```
haplotype_combined.csv
```

in your chosen output directory.

------

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

