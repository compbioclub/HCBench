##  🧬 HCBench Parsers 

Welcome to the **HCBench Parser Suite** — a unified framework for standardizing single-cell copy number alteration (CNA) outputs from multiple tools.

Each parser converts heterogeneous output formats from different algorithms (CHISEL, Alleloscope, CNRein, SEACON, and SIGNALS) into a single consistent format suitable for downstream benchmarking and visualization.

------

### 🎯 Common Standardized Outputs

While each tool requires different input files and formats, every parser in this suite is designed to output a canonical directory structure. Depending on the specific methods called, a fully processed output directory will typically contain:

**CNA Matrices**: `haplotype_combined.csv`, or split `minor.csv`, `major.csv`, and `minor_major.csv` matrices (regions × cells).

**Cluster Mapping**: `clusters.csv` containing standardized `cell_id` and `clone_id` columns.

**Bin Counts / RDR**: `bin_counts.csv` or `bin_rdr.csv` formatted as a wide matrix of regions by cells.

**Sparse VAF Matrices**: A `VAF/` directory containing standard Matrix Market (`.mtx`) files for Allelic Depth (AD) and Read Depth (DP).



## 📖 Parsers Overview

| **Parser**                                                   | **Main Input**                  | **Additional Inputs**                             | **Format Types**            | **Supported Outputs**                        |
| ------------------------------------------------------------ | ------------------------------- | ------------------------------------------------- | --------------------------- | -------------------------------------------- |
| **[CHISEL](chisel.md)**      | `calls.tsv`                     | `mapping.tsv`, VAF table                          | Tab-delimited text          | CNA Matrix, Clusters, Bin Counts, Sparse VAF |
| **[Alleloscope](alleloscope.md)** | `.rds` files                    | `clusters.csv`, raw counts TSV, cellSNP directory | R serialized data, TSV, CSV | CNA Matrix, Clusters, Bin Counts, Sparse VAF |
| **[CNRein](cnrein.md)**      | `CNReinPrediction.csv`          | `.npz` arrays, split VCF files                    | CSV, NPZ, VCF               | CNA Matrix, Bin RDR, Sparse VAF              |
| **[SEACON](seacon.md)**      | `calls.tsv`                     | `counts.tsv`, `vaf.tsv`                           | Tab-delimited text          | Split CNA Matrices, Bin Counts, Sparse VAF   |
| **[SIGNALS](signals.md)**                                                  | `hscn.rds` (exported to `.tsv`) | cluster file, bin counts, VAF table               | Tab-delimited text / CSV    | CNA Matrix, Clusters, Bin Counts, Sparse VAF |

All parsers generate outputs in the same canonical structure, making it possible to directly compare results across tools in the **HCBench** benchmarking pipeline.