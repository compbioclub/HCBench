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

## Manual CNA annotation

CNA annotation can be performed using the `annotate_segments()` method provided by each tool parser. This method annotates standardized CNA profiles with segment information, event sizes, and CNA types.

### Annotate CNA profiles

After generating standardized CNA profiles with `parser.run()`, use `annotate_segments()` to annotate the desired output file:

```python
parser.run()

annotated = parser.annotate_segments(
    cna_df_path="/results/chisel/cell_level/haplotype_combined.csv",
    out_file="CHISEL_annotated.csv",
    detect_cnv_type=True,
    # threshold=0.8,
    # small_size=3_000_000,
    # mid_size=10_000_000,
)
```

By default, the method uses a threshold of `0.8` and classifies CNA events as focal (less than 3 Mb), medium (3–10 Mb), or broad (greater than 10 Mb).

The `genome` parameter specifies the reference genome (default: the parser's `ref_genome`), while `chromosomes` allows annotation of selected chromosomes (default: all chromosomes).

The annotated results are saved to `<parser.output_path>/<out_file>` (default: `annotated_segments.csv`) and returned as a DataFrame. When specific chromosomes are provided, additional chromosome-specific files are generated with suffixes such as `_chr1.csv`.

For hcCNA analysis, use haplotype-specific CN profiles rather than unphased minor/major CN values.

### Annotate ground truth CNA profiles

Ground truth (GT) CNA profiles can be annotated using `GTParser.annotate_segments()`:

```python
gt_parser.annotate_segments(
    "/results/gt/haplotype_combined.csv",
    threshold=0.9,
    detect_cnv_type=True,
    out_file="ground_truth_classified.csv",
)
```

The default threshold for GT annotation is `0.9`, and `detect_cnv_type` is disabled by default. Set `detect_cnv_type=True` to include CNA event types in the output, particularly when preparing annotated files for complex hcCNA comparisons.

The annotation identifies CNA segments, assigns event sizes and types, and includes qualifying whole-chromosome duplication (WCD), whole-chromosome loss (WCL), and whole-genome duplication (WGD) events.

For evaluating complex hcCNA events using annotated files, see [RealBench](../realbench/index.md).