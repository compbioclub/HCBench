# Alleloscope Parser

This module provides `AlleloscopeParser`, specialized for parsing **Alleloscope** outputs and exporting standardized matrices and helper files used by hcbench workflows.

Key features:

- Read Alleloscope **RDS** outputs (genotypes and segmentation table) and convert them into a unified haplotype level CNA matrix
- Optional region splitting by a user defined `bin_size`
- Parse Alleloscope cluster assignments into a standardized `clusters.csv`
- Export bin-level count matrices as `bin_counts.csv`
- Convert a cellSNP-like VAF long table into sparse matrix outputs

## 🚀 Quick Start

### 🧬 1. Parse the CNA Matrix from RDS

```python
from hcbench.parsers.alleloscope import AlleloscopeParser

allelo_output = "/output/alleloscope/"
genotypes_rds = "/demo_output/alleloscope/output/rds/genotypes.rds"
seg_table_rds = "/demo_output/alleloscope/output/rds/seg_table.rds"

alleloscope_parser = AlleloscopeParser(
    output_path=allelo_output,
    genotypes_rds_path=genotypes_rds,
    seg_table_rds_path=seg_table_rds,
    barcode_path=None,
    bin_size=100000,          # set an integer to split regions, e.g. 100000
    start_offset=1,         # shift start coordinate by +1 if needed
)

alleloscope_parser.run()
```

After running, the parser will read the two RDS files and save results to the output directory, typically containing files, for example:

```
/output/alleloscope/
  haplotype_combined.csv
  haplotype_1.csv
  haplotype_2.csv
  minor.csv
  major.csv
  minor_major.csv
```

### 🧬 2. Parse the Cluster File

```python
cluster_file = "/demo_output/alleloscope/output/clusters.csv"
alleloscope_parser.get_cluster(cluster_file)
```

After running this command, the following standardized file will be created:

```
/output/alleloscope/
  clusters.csv
```

with content:

```
cell_id,clone_id
cellA,1
cellB,1
cellC,2
```

### 🧬 3. Parse the bin counts matrix

```
counts_file = "/demo_output/alleloscope/input/counts/raw_counts.tsv"
alleloscope_parser.get_bin_counts(counts_file)
```

After running this command, the following standardized file will be created:

```
/output/alleloscope/
  bin_counts.csv
```

### 🧬 4. Parse the VAF sparse matrices

This method assumes `vaf_file_path` is a directory containing cellSNP style outputs such as:

- `cellSNP.samples.tsv`
- `cellSNP.tag.AD.mtx`
- `cellSNP.tag.DP.mtx`
- `cellSNP.base.vcf.gz`

```
alleloscope_parser.get_VAF_matrix(
    vaf_file_path="/demo_output/alleloscope/cellSNP_out",
    min_dp=3,
    min_cells=10,
    prefix="cellSNP",
)
```

Output:

```
/output/alleloscope/
  VAF/
    cellSNP.VAF.filtered.mtx
    ...
```



## ⚙️ Initialization

```
AlleloscopeParser(
    output_path: str,
    genotypes_rds_path: str,
    seg_table_rds_path: str,
    barcode_path: str | None = None,
    bin_size: int | None = None,
    add_chr_prefix: bool = True,
    start_offset: int = 1,
    **kwargs
)
```

**`output_path`** : base output directory.

**`genotypes_rds_path`** : Path to the Alleloscope genotypes RDS file.

**`seg_table_rds_path`**: Path to the Alleloscope segmentation table RDS file. It is used to map region identifiers into standardized region labels like:

- `chr<CHR>:<START+start_offset>-<END>` when `add_chr_prefix=True`
- `<CHR>:<START+start_offset>-<END>` when `add_chr_prefix=False`

**`barcode_path`** (optional) :Path to a barcode mapping file used to remap cell IDs

**`bin_size`** (optional) :If set, regions are split into fixed sized. 

**`add_chr_prefix`** (default `True`) :Controls whether to prepend `"chr"` when forming the region label from the segmentation table.

**`start_offset`** (default `1`) :Adds an offset to the start coordinate when forming region labels.

## 🧠 Core Methods









## 🧬 Generating a Cluster File in R

If you wish to generate the Alleloscope clustering file from your processed object in R,
 you can use the following commands:

```
linplot = Lineage_plot(Obj_filtered = Obj_filtered, nSNP = 2000, nclust = 10)
write.csv(linplot, file = "cluster.csv", row.names = TRUE)
```

This will create a `cluster.csv` file that can be used as input to `AlleloscopeParser.get_cluster()` in Python.

------
