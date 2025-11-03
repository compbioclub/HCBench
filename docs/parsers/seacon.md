
# SEACON Parser

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
