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
the parser automatically **increments the `START` value by 1** — converting 0-based to 1-based coordinates, as per standard genomic conventions.

---

## 📤 Output Files

After parsing, results are saved to the `chisel_output/` directory, containing the following six files:

```
haplotype_combined.csv
haplotype_1.csv
haplotype_2.csv
minor.csv
major.csv
minor_major.csv
```



---

## 🚀 Example Usage

```python
from hcbench.parsers.chisel import ChiselParser

chisel_input = "/mnt/cbc_adam/public/workspace/HCDSIM/demo_output/chisel/calls/calls.tsv"
chisel_output = "/home/jianganna/workspace/HCBench_project/HCBench/output/chisel/"

chisel_parser = ChiselParser(chisel_input, chisel_output)
chisel_parser.run()
```

After running, the parser will read `calls.tsv` and automatically generate the standardized output files listed above.