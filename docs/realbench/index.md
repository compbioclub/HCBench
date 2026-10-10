# HCBench-real

HCBench-real evaluates copy number alteration (CNA) results from real datasets by comparing predictions from different callers. It supports five analyses: CNA detection, clone detection, evolution tracking, accumulated complex hcCNA detection, and CN phasing.

As ground truth is generally unavailable for real datasets, most metrics describe **agreement between callers**, rather than accuracy against a known truth. The read-depth (RDR) and variant allele frequency (VAF) analyses instead compare inferred CN profiles with observed sequencing signals.

The examples below compare Alleloscope and CHISEL. Replace the input paths and `tools` list with your own results. `acCNA` refers to allele- and cell/clone-specific CNAs, while `hcCNA` refers to haplotype- and cell/clone-specific CNAs.

## Setup

Run the examples from the repository root, or set `example_dir` to the location of your `exp/` directory. Input files remain separate from analysis outputs. Both `exp/Alleloscope/` and `exp/CHISEL/` contain `bin_counts.csv`, `clusters.csv`, `haplotype_1.csv`, `haplotype_2.csv`, `haplotype_combined.csv`, `major.csv`, `minor.csv`, and `minor_major.csv`. All input lists below follow the same tool order.

```python
from pathlib import Path
from hcbench.realbench import RealBench

example_dir = Path("exp")
output_dir = "output/realbench_examples"
realbench_runner = RealBench(output_dir=output_dir)

tools = ["Alleloscope", "CHISEL"]
```

------

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
)
```

By default, `complex_hccnas()` annotates the input CN profiles and compares complex hcCNA calls between tools. Use `genome="hg38"` (default) or `genome="hg19"` to match the reference genome of your data. The default annotation threshold is `0.8`.

To save the annotated events for each tool, set `save_annotated_df=True`. The files are written to `{output_dir}/{outprefix}annotated/` as `{tool}_annotated.csv`. You can also use `chromosomes=["chr1", "chr2"]` to limit the analysis to selected chromosomes; both `"chr1"` and `"1"` formats are accepted. By default, all chromosomes are included.

#### Use existing annotated CNA files

If your CN profiles have already been annotated, provide the annotation tables directly to avoid repeating this step:

To reuse previously annotated events, first generate them with `save_annotated_df=True` (or through the Parser). Supply the resulting file paths in `annotated_files`, keyed by caller. For the example above, annotations are saved under `{output_dir}/complex_hccnas_0.8annotated/`:

```python
realbench_runner.complex_hccnas(
    annotated_files=annotated_files,
    outprefix="complex_hccnas_reused_",
)
```

Each tool requires one combined CSV file with the columns `cell`, `Chrom`, `Start`, `End`, `value`, `Segment`, `size`, and `type`. Coordinates must be positive, 1-based, and inclusive. The files should contain the complete independently annotated events, rather than only focal events or events derived from pairwise-aligned data. You can also pass a list of file paths through `annotated_files`, together with matching `tool_names`. Do not use `tool_cna_files` and `annotated_files` in the same call.

The comparison aligns event intervals using shared genomic boundaries while retaining each tool's original event types, sizes, and CN values. Regions without an annotated event are not automatically classified as neutral. The resulting `overlap.csv` and `metrics.csv` report pairwise agreement across shared cell/region entries and event categories; these metrics are not base-pair-weighted measures of event overlap.

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

#### Prepare gene and pathway annotations

This analysis requires annotated hcCNA events, a gene annotation (GTF) file, and pathway gene sets. HCBench includes three pathway libraries: `MSigDB_Hallmark_2020`, `Reactome_2022`, and `KEGG_2021_Human`.

For gene annotation, HCBench uses [GENCODE v46 (hg38)](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz) or [GENCODE v46lift37 (hg19)](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/gencode.v46lift37.annotation.gtf.gz). If `gtf_path` is not specified, the appropriate GTF file is downloaded automatically when needed and cached in `~/.cache/hcbench/resources/genes/<genome>/`. Alternatively, provide the path to a local GTF file using `gtf_path`. Make sure `genome` matches the genomic coordinates in your CNA files.

By default, all three bundled pathway libraries are available. Use `libraries` to select specific libraries, `gene_sets` to provide custom GMT files, or `resource_dir` to use a custom resource directory. When using the bundled libraries, you can leave `resource_dir` unset.

#### Run focal genes and pathways analysis

Generate annotated hcCNA events before running this analysis. Set `annotated_inputs` to records containing `sample`, `tool`, `bin_size`, and `files` for each caller and resolution. To run offline, set `gtf_path` to a locally available GTF file:

```python
# Focal genes and pathways
annotated_inputs=[
            {
                "sample": "sample1",
                "tool": "Alleloscope",
                "bin_size": "5Mb",
                "files": f"{output_dir}/complex_hccnas_0.8annotated/Alleloscope_annotated.csv"
            },
            {
                "sample": "sample1",
                "tool": "CHISEL",
                "bin_size": "5Mb",
                "files": [
                    f"{output_dir}/complex_hccnas_0.8annotated/CHISEL_annotated.csv",
                ]
            }
        ]
realbench_runner.focal_genes_and_pathways(
    annotated_inputs=annotated_inputs,
    genome="hg38",
    gtf_path=None,  # optionally replace None with a local GTF path
    # outprefix="focal_genes_and_pathways"
)
```

Each entry in `annotated_inputs` identifies a sample, a caller, a bin size, and its annotated CNA file(s). You can provide either a combined CSV covering multiple chromosomes or a list of chromosome-specific CSV files. Multiple entries allow comparisons across callers and resolutions.

The analysis reports genes overlapping focal hcCNA events, over-representation analysis (ORA) results, gene- and pathway-level Jaccard similarities, and summary statistics. See the output files for results from each sample, caller, and resolution.

------

## CN Phasing

Compare genome-wide CN phasing between callers at the cell level.

### Global CN phasing

Calculate mismatch error between inferred CN phases in shared cells. With `mode="heterozygous-only"`, comparisons are restricted to heterozygous loci informative in both callers; lower mismatch error indicates greater agreement.

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
    mode = "heterozygous-only",
    # outprefix = "global_cn_phasing_heterozygous-only"
)
```
