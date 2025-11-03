


# CNRein Parser

## 📂 Input Files

The output directory of **CNRein** typically includes a single main file that records per-cell haplotype-specific copy number information:

```
demo_output/cnrein/finalPrediction
└── CNReinPrediction.csv
```

- `CNReinPrediction.csv` — main output file containing haplotype-level copy number information for each cell across all genomic regions.

An example of `CNReinPrediction.csv`:

```
Cell barcode,Chromosome,Start,End,Haplotype 1,Haplotype 2
TN1_3_S1_C85,1,800001,20100000,2,0
TN1_3_S1_C85,1,20100001,29200000,3,1
TN1_3_S1_C85,1,29200001,33100000,3,1
TN1_3_S1_C85,1,33100001,34900000,3,1
```

The **required columns** are:

```
Chromosome, Start, End, Cell barcode, Haplotype 1, Haplotype 2
```

The parser combines the two haplotype columns (`Haplotype 1` and `Haplotype 2`) into a single unified field called `HAP_CN`,
 formatted as `"hap1|hap2"`, for example `"2|0"`.

------

## 📤 Output Files

After parsing, the following file is generated in the specified output directory:

```
haplotype_combined.csv
```

- `haplotype_combined.csv` — the standardized CNA matrix (regions × cells).
   Each entry represents the combined haplotype state of a cell in `"Haplotype 1|Haplotype 2"` format.

If the `split_haplotype=True` option is used in the parser base class,
 additional per-haplotype and major/minor copy number matrices will also be generated automatically:

```
haplotype_1.csv
haplotype_2.csv
minor.csv
major.csv
minor_major.csv
```

------

## ⚙️ Key Parameters

| Parameter        | Description                                           | Default          |
| ---------------- | ----------------------------------------------------- | ---------------- |
| `chrom_col`      | Column name for chromosome.                           | `"Chromosome"`   |
| `start_col`      | Column name for region start.                         | `"Start"`        |
| `end_col`        | Column name for region end.                           | `"End"`          |
| `cell_col`       | Column name for cell barcode.                         | `"Cell barcode"` |
| `value_col`      | Column name for the combined haplotype CN value.      | `"HAP_CN"`       |
| `start_plus_one` | Whether to convert 0-based to 1-based coordinates.    | `False`          |
| `add_chr_prefix` | Whether to ensure `"chr"` prefix on chromosome names. | `True`           |

------

## 🚀 Example Usage

```
from hcbench.parsers.cnrein import CNReinParser

cnrein_input = "/demo_output/cnrein/finalPrediction/CNReinPrediction.csv"
cnrein_output = "/output/cnrein/"

cnrein_parser = CNReinParser(
    input_path=cnrein_input,
    output_path=cnrein_output
)

cnrein_parser.run()
```

This will read `CNReinPrediction.csv`, merge `Haplotype 1` and `Haplotype 2` into the `HAP_CN` column,
 and save a standardized file named:

```
haplotype_combined.csv
haplotype_1.csv       # if split_haplotype=True
haplotype_2.csv       # if split_haplotype=True
minor.csv             # if split_haplotype=True
major.csv             # if split_haplotype=True
minor_major.csv       # if split_haplotype=True
```

in your output directory.

------

