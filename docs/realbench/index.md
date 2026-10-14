# HCBench-real

HCBench-real evaluates copy number alteration (CNA) results from real datasets by comparing predictions from different callers. It supports five analyses: CNA detection, clone detection, evolution tracking, accumulated complex hcCNA detection, and CN phasing.

As ground truth is generally unavailable for real datasets, most metrics describe **agreement between callers**, rather than accuracy against a known truth. The read-depth (RDR) and variant allele frequency (VAF) analyses instead compare inferred CN profiles with observed sequencing signals.

The examples below compare Alleloscope and CHISEL. Replace the input paths and `tools` list with your own results. `acCNA` refers to allele- and cell/clone-specific CNAs, while `hcCNA` refers to haplotype- and cell/clone-specific CNAs.

## Prepare inputs and initialize

`bin_size` and `genome` are required when creating a CNA parser, `GTBench`, or `RealBench`:

- `bin_size`: a positive integer in base pairs, such as `100000` for 100 kb. Use the resolution of your CNA data.
- `genome`: `"hg19"` or `"hg38"`, matching your input coordinates.

Methods use this configuration. Initialization prints:

```text
Bin size: 100000 bp
Genome: hg38
```

With the default event-size thresholds, analysis is available for:

| CNA `bin_size` | Event categories |
| --- | --- |
| `<= 3,000,000` bp | focal, medium, broad |
| `> 3,000,000` and `<= 10,000,000` bp | medium, broad |
| `> 10,000,000` bp | broad |

Event length is `End - Start + 1`: focal is below 3 Mb, medium is from 3 Mb through 10 Mb, and broad is above 10 Mb. Whole-chromosome and whole-genome events are reported as broad. Categories unavailable at your resolution have no annotation or category-metric rows; a missing comparison is reported as NA rather than evidence of no events.

The category limits apply to complex-event annotation and mirrored-event analysis. If `bin_size > 3,000,000`, focal gene/pathway analysis prints a skip message, returns `None`, and writes no new enrichment results.

Run this setup once before the method examples. Standardized CNA matrices have a first `region` column and remaining cell columns. Clone files have `cell_id` and `clone_id`; all callers should describe the same sample and genome assembly.

Run the examples from the repository root, or set `example_dir` to the location of your `exp/` directory. Input files remain separate from analysis outputs. Both `exp/Alleloscope/` and `exp/CHISEL/` contain `bin_counts.csv`, `clusters.csv`, `haplotype_1.csv`, `haplotype_2.csv`, `haplotype_combined.csv`, `major.csv`, `minor.csv`, and `minor_major.csv`. All input lists below follow the same tool order.

```python
from pathlib import Path
from hcbench.realbench import RealBench

example_dir = Path("exp")
output_dir = "output/realbench_examples"
realbench_runner = RealBench(output_dir=output_dir, bin_size=100000, genome="hg38")

tools = ["Alleloscope", "CHISEL"]
```

------

Outputs are saved under `output_dir`; it is created during initialization. `cna_values()` writes `{outfile_prefix}_rmse.csv`, `{outfile_prefix}_scc.csv`, and `{outfile_prefix}_acc.csv` as caller-by-caller matrices and returns `None`. Clone assignment comparisons return ARI/AMI DataFrames and save CSV matrices. Sections below describe inputs needed beyond CNA profiles.

## CNA Detection

Compare allele-specific and haplotype-specific CN values and states between callers. The RDR and VAF analyses additionally measure agreement with observed sequencing data.

### acCNA values

Compare allele-specific CN values at shared cells and genomic loci using root mean square error (RMSE), Spearman correlation coefficient (SCC), and accuracy (ACC).

```python
# acCNA values
realbench_runner.cna_values(
        tool_cna_files=[
            f"{example_dir}/Alleloscope/minor_major.csv",
            f"{example_dir}/CHISEL/minor_major.csv",
        ],
        tool_names=tools,
        haplotype = "combined",
        outfile_prefix= "cna_values_accna"
    )
```

### acCNA states

Compare allele-specific CN loss, gain, and neutral states at shared cells and loci using ACC and Cohen's kappa.

```python
# acCNA states
realbench_runner.cna_states(
    tool_hap1_cna_files=[
        f"{example_dir}/Alleloscope/minor.csv",
        f"{example_dir}/CHISEL/minor.csv",
    ],
    tool_hap2_cna_files=[
        f"{example_dir}/Alleloscope/major.csv",
        f"{example_dir}/CHISEL/major.csv",
    ],
    tool_names=tools,
    outfile="cna_states_accna.csv",
)
```

### hcCNA values

Compare haplotype-specific CN values at shared cells and loci using RMSE, SCC, and ACC. Haplotype order is retained when comparing phased CN values.

```python
# hcCNA values
realbench_runner.cna_values(
        tool_cna_files=[
            f"{example_dir}/Alleloscope/haplotype_combined.csv",
            f"{example_dir}/CHISEL/haplotype_combined.csv",
        ],
        tool_names=tools,
        haplotype = "combined",
        outfile_prefix= "cna_values_hccna"
    )
```

### hcCNA states

Compare haplotype-specific CN loss, gain, and neutral states at shared cells and loci using ACC and Cohen's kappa.

```python
# hcCNA states
realbench_runner.cna_states(
    tool_hap1_cna_files=[
        f"{example_dir}/Alleloscope/haplotype_1.csv",
        f"{example_dir}/CHISEL/haplotype_1.csv",
    ],
    tool_hap2_cna_files=[
        f"{example_dir}/Alleloscope/haplotype_2.csv",
        f"{example_dir}/CHISEL/haplotype_2.csv",
    ],
    tool_names=tools,
    outfile="cna_states_hccna.csv",
)
```

### Total CN vs. RDR

Measure the agreement between inferred total CN and observed read-depth ratios (RDR) using L1 error (mean absolute error). Lower values indicate closer agreement. Paired t-tests compare within-cell L1 errors on cells shared by each pair of callers.

```python
# Total CN vs. RDR
realbench_runner.total_cn_vs_rdr(
        bin_count_files=[
            f"{example_dir}/Alleloscope/bin_counts.csv",
            f"{example_dir}/CHISEL/bin_counts.csv",
        ],
        tool_cna_files=[
            f"{example_dir}/Alleloscope/haplotype_combined.csv",
            f"{example_dir}/CHISEL/haplotype_combined.csv",
        ],
        tool_names=tools,
)
```

### Major/minor CN vs. VAF

Assess the fit of inferred major/minor CN distributions to observed variant allele frequencies (VAFs) with a binomial-model log-likelihood ratio (LLR).

This example dataset does not include SNV matrices. To run this analysis, provide `snv_paths` containing `cellSNP.AD.filtered.mtx` and `cellSNP.DP.filtered.mtx` for each caller, along with `cell_lists` (matrix columns) and `variant_pos_dfs` (matrix rows). Keep all lists in `tools` order:

```python
# Major/minor CN vs. VAF
realbench_runner.major_minor_cn_vs_vaf(
    tool_hap1_files=[
        f"{example_dir}/Alleloscope/haplotype_1.csv",
        f"{example_dir}/CHISEL/haplotype_1.csv",
    ],
    tool_hap2_files=[
        f"{example_dir}/Alleloscope/haplotype_2.csv",
        f"{example_dir}/CHISEL/haplotype_2.csv",
    ],
    tool_names=tools,
    snv_paths=snv_paths,
    cell_lists  = cell_lists,
    variant_pos_dfs = variant_pos_dfs,
)
```

------

## Clone Detection

Compare clone assignments, clone sizes, and unique CN profiles between callers.

### Global tumor clones

Compare cell-to-clone assignments using adjusted mutual information (AMI) and adjusted Rand index (ARI).

```python
# Global tumor clones
realbench_runner.global_tumor_clones(
        tool_cluster_files=[
            f"{example_dir}/Alleloscope/clusters.csv",
            f"{example_dir}/CHISEL/clusters.csv",
        ],
        tool_names=tools,
)
```

### Rare and dominant clones

Use clone size deviation (CSD) to compare the sizes of clones assigned to shared cells, including rare and dominant clones.

```python
# Rare and dominant clones
realbench_runner.rare_and_dominant_clones(
        tool_cluster_files=[
            f"{example_dir}/Alleloscope/clusters.csv",
            f"{example_dir}/CHISEL/clusters.csv",
        ],
        tool_names=tools,
)
```

### Unique acCNA profiles

Compare the number and size distribution of unique allele-specific CN profiles using unique profile count (UPC) and unique profile size deviation (UPSD).

```python
# Unique acCNA profiles
realbench_runner.unique_cna_profile_counts(
        tool_cna_files=[
            f"{example_dir}/Alleloscope/minor_major.csv",
            f"{example_dir}/CHISEL/minor_major.csv",
        ],
        tool_names=tools,
    outfile="unique_cna_profile_counts_accna.csv",
    )

realbench_runner.unique_cna_profile_sizes(
        tool_cna_files=[
            f"{example_dir}/Alleloscope/minor_major.csv",
            f"{example_dir}/CHISEL/minor_major.csv",
        ],
        tool_names=tools,
    outfile="unique_cna_profile_sizes_accna.csv",
    )
```

### Unique hcCNA profiles

Compare the number and size distribution of unique haplotype-specific CN profiles using UPC and UPSD.

```python
# Unique hcCNA profiles
realbench_runner.unique_cna_profile_counts(
        tool_cna_files=[
            f"{example_dir}/Alleloscope/haplotype_combined.csv",
            f"{example_dir}/CHISEL/haplotype_combined.csv",
        ],
        tool_names=tools,
    outfile="unique_cna_profile_counts_hccna.csv",
    )

realbench_runner.unique_cna_profile_sizes(
        tool_cna_files=[
            f"{example_dir}/Alleloscope/haplotype_combined.csv",
            f"{example_dir}/CHISEL/haplotype_combined.csv",
        ],
        tool_names=tools,
    outfile="unique_cna_profile_sizes_hccna.csv",
    )
```

------

## Evolution Tracking

Evaluate phylogenetic parsimony scores (PS) derived from caller-specific CN profiles.

### acCNA phylogeny

Calculate parsimony scores (PS) from allele-specific CN profiles to compare inferred evolutionary complexity between callers.

```python
# acCNA phylogeny
realbench_runner.cna_phylogeny(
        tool_cna_files=[
            f"{example_dir}/Alleloscope/minor_major.csv",
            f"{example_dir}/CHISEL/minor_major.csv",
        ],
        tool_names=tools,
    outfile="cna_phylogeny_accna.csv",
    )
```

### hcCNA phylogeny

Calculate parsimony scores (PS) from haplotype-specific CN profiles to compare inferred evolutionary complexity between callers.

```python
# hcCNA phylogeny
realbench_runner.cna_phylogeny(
        tool_cna_files=[
            f"{example_dir}/Alleloscope/haplotype_combined.csv",
            f"{example_dir}/CHISEL/haplotype_combined.csv",
        ],
        tool_names=tools,
    outfile="cna_phylogeny_hccna.csv",
    )
```

------

## Accumulated Complex hcCNA Detection

Compare accumulated complex hcCNA and mirrored-subclonal-CNA events at the cell level. For focal complex hcCNAs, compare affected genes and pathways at the sample level.

### Complex hcCNAs

Compare accumulated complex hcCNA events between callers at the cell level. Report Jaccard similarity (JS) for event overlap and RMSE, ACC, and SCC for CN-value agreement.

```python
# Complex hcCNAs
realbench_runner.complex_hccnas(
    tool_cna_files=[
        f"{example_dir}/Alleloscope/haplotype_combined.csv",
        f"{example_dir}/CHISEL/haplotype_combined.csv",
    ],
    tool_names=tools,
    threshold=0.8,
    outprefix="complex_hccnas_0.8",
    save_annotated_df=True,
)
```

`complex_hccnas()` writes `complex_hccnas_0.8overlap.csv` and `complex_hccnas_0.8metrics.csv` in `output_dir` for this example, and returns `None`. The annotation threshold defaults to `0.8`. Missing category comparisons are NA; categories unavailable at the configured resolution are omitted.

To save the annotated events for each tool, set `save_annotated_df=True`. The files are written to `{output_dir}/{outprefix}annotated/` as `{tool}_annotated.csv`. You can also use `chromosomes=["chr1", "chr2"]` to limit the analysis to selected chromosomes; both `"chr1"` and `"1"` formats are accepted. By default, all chromosomes are included.

#### Use existing annotated CNA files

Use the annotation CSVs saved by the preceding call or produced by a parser. Supply one file per caller:

```python
annotated_files = {
    tool: f"{output_dir}/complex_hccnas_0.8annotated/{tool}_annotated.csv"
    for tool in tools
}
realbench_runner.complex_hccnas(
    annotated_files=annotated_files,
    outprefix="complex_hccnas_reused_",
)
```

Each tool requires one combined CSV file with the columns `cell`, `Chrom`, `Start`, `End`, `value`, `Segment`, `size`, and `type`. Coordinates must be positive, 1-based, and inclusive. The files should contain the complete independently annotated events, rather than only focal events or events derived from pairwise-aligned data. You can also pass a list of file paths through `annotated_files`, together with matching `tool_names`. Do not use `tool_cna_files` and `annotated_files` in the same call.

The comparison aligns event intervals using shared genomic boundaries using each tool's event types, sizes, and CN values. Regions without an annotated event are not automatically classified as neutral. The resulting `overlap.csv` and `metrics.csv` report pairwise agreement across shared cell/region entries and event categories; these metrics are not base-pair-weighted measures of event overlap.

For instructions on generating annotated CNA tables, see [Parser annotation](../parsers/index.md).

### Mirrored-subclonal-CNAs

Compare accumulated mirrored-subclonal-CNA events between callers at the cell level using JS, RMSE, and ACC.

```python
# Mirrored-subclonal-CNAs
realbench_runner.mirrored_subclonal_cnas(
        tool_cna_files=[
            f"{example_dir}/Alleloscope/haplotype_combined.csv",
            f"{example_dir}/CHISEL/haplotype_combined.csv",
        ],
        tool_names=tools,
        use_segemnt = True,
        save_tmp = True
    )
```


### Focal genes and pathways

Identify genes overlapping focal complex hcCNA events, run pathway over-representation analysis (ORA), and compare affected gene and pathway sets between callers and resolutions using Jaccard similarity (JS).

This analysis requires the `RealBench` initialization `bin_size` to be at most the focal threshold, currently 3,000,000 bp. At exactly 3 Mb it is allowed; above 3 Mb the method prints the configured bin size and threshold, skips analysis, returns `None`, and creates no new enrichment results.

Install the `enrichment` extra and use the annotation files saved by the complex-event example above. Their resolution labels below match the 100 kb initialization. HCBench downloads and caches a matching GTF when needed; to run offline, provide a local `gtf_path`.

```python
annotated_inputs = [
    {
        "sample": "sample1",
        "tool": tool,
        "bin_size": "100kb",
        "files": f"{output_dir}/complex_hccnas_0.8annotated/{tool}_annotated.csv",
    }
    for tool in tools
]
result = realbench_runner.focal_genes_and_pathways(
    annotated_inputs,
    workers=1,
)
if result is not None:
    print(result["output_dir"])
```

Results are saved under `output_dir/focal_genes_and_pathways/`, including `group_summary.csv`, `focal_pathways_significant.csv`, `consistency_jaccard.csv`, and `enrichment_status.csv`. Read the status file before interpreting an absent pathway result. See the [focal gene and pathway tutorial](focal_genes_and_pathways.md) for input requirements, local references, background choices, and return values.

------

## CN Phasing

Compare genome-wide CN phasing between callers at the cell level.

### Global CN phasing

Calculate mismatch error between inferred CN phases in shared cells. With `mode="Shared-heterozygous"`, comparisons are restricted to heterozygous loci informative in both callers; lower mismatch error indicates greater agreement.

```python
# Global CN phasing
realbench_runner.global_cn_phasing(
    tool_hap1_cna_files=[
        f"{example_dir}/Alleloscope/haplotype_1.csv",
        f"{example_dir}/CHISEL/haplotype_1.csv",
    ],
    tool_hap2_cna_files=[
        f"{example_dir}/Alleloscope/haplotype_2.csv",
        f"{example_dir}/CHISEL/haplotype_2.csv",
    ],
    tool_names=tools,
    mode = "Shared-heterozygous",
    # outprefix = "global_cn_phasing_Shared-heterozygous"
)
```
