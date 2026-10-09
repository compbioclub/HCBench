# HCBench-real

`HCBench-real` is an integrated module for evaluating copy number alteration (CNA) results in **real-world datasets**. Since real data often lacks a ground truth, this module focuses on pairwise cross-tool consistency and agreement with observed sequencing signals. The evaluation is organized around five evaluation questions: CNA detection, clone detection, evolution tracking, accumulated complex hcCNA detection, and CN phasing.

## Setup

Run the examples from the repository root, or set `example_dir` to the location of your `exp/` directory. Input files remain separate from analysis outputs. Both `exp/Alleloscope/` and `exp/CHISEL/` contain `bin_counts.csv`, `clusters.csv`, `haplotype_1.csv`, `haplotype_2.csv`, `haplotype_combined.csv`, `major.csv`, `minor.csv`, and `minor_major.csv`. All input lists below follow the same tool order.

```python
from pathlib import Path
from hcbench.realbench.realbench import RealBench

example_dir = Path("exp")
output_dir = "output/realbench_examples"
realbench_runner = RealBench(output_dir=output_dir)

tools = ["Alleloscope", "CHISEL"]
```

------

## CNA Detection

This section assesses pairwise consistency of allele-specific and haplotype-specific CN calls and evaluates whether predicted copy numbers agree with observed read depth and allele-frequency signals.

### acCNA values

**Are allele-, cell-, and locus-specific CN values consistent?**

Quantifies the direct numerical agreement of copy number values across all shared cells and loci using RMSE, SCC, and ACC.

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

**Are allele-, cell-, and locus-specific CN loss, gain, and neutral states consistent?**

Uses a binary classification framework to measure concordance. By treating one caller as a reference, metrics like ACC and the Kappa score are used to determine how often the two tools agree on the presence of gains, losses, or neutral states.

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

**Are haplotype-, cell-, and locus-specific CN values consistent?**

Measures numerical consistency for phased copy numbers (e.g., matching 1|2 vs 1|2) using RMSE, SCC, and the exact match rate (ACC).

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

**Are haplotype-, cell-, and locus-specific CN loss, gain, and neutral states consistent?**

Evaluates the reliability of phased CN calling by measuring the concordance (ACC and Kappa score) between callers when categorizing segments into gain, loss, or neutral states.

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

**Are total CNs consistent with observed read depths?** Evaluates the goodness of fit between the inferred total copy number and the observed read-depth signal using the Read Depth L1 Error. Statistically compares callers using paired t-tests; a lower L1 error indicates that the predicted states more accurately reflect the raw data counts.

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

**Are major and minor CN distributions consistent with observed VAFs?** Measures how well the predicted major and minor alleles align with truncal SNV frequencies (or BAFs) using log-likelihood under a binomial model. A higher log-likelihood ratio (LLR) indicates a superior fit to the observed molecular evidence.

`exp/` contains the CNA inputs but no VAF/SNV matrices. Before using this API, provide `snv_paths` for Alleloscope and CHISEL, each containing `cellSNP.AD.filtered.mtx` and `cellSNP.DP.filtered.mtx`, plus `cell_lists` matching their matrix columns and `variant_pos_dfs` matching their variant rows. All three lists must follow `tools` order. No VAF paths are assumed to exist in `exp/`; the following is a conditional call template:

```python
# Major/minor CN vs. VAF
# realbench_runner.major_minor_cn_vs_vaf(
#     tool_hap1_files=[
#         f"{example_dir}/Alleloscope/haplotype_1.csv",
#         f"{example_dir}/CHISEL/haplotype_1.csv",
#     ],
#     tool_hap2_files=[
#         f"{example_dir}/Alleloscope/haplotype_2.csv",
#         f"{example_dir}/CHISEL/haplotype_2.csv",
#     ],
#     tool_names=tools,
#     snv_paths=snv_paths,
#     cell_lists  = cell_lists,
#     variant_pos_dfs = variant_pos_dfs,
# )
```

------

## Clone Detection

This section evaluates pairwise consistency of inferred clone architectures and the diversity of unique CN profiles between callers.

### Global tumor clones

**Are global tumor clones consistent?**

Evaluates the agreement of clustering labels between two callers using AMI and ARI. These metrics quantify how similar the overall cell groupings are across the entire population.

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

**Are rare and dominant clones consistent?**

Assesses consistency stratified by predicted clone size. Using Clone Size Deviation (CSD), it examines whether a cell assigned to a cluster of size $n$ by one caller is assigned to a cluster of similar size by the other, identifying if callers disagree on the granularity of specific subclones.

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

**Are unique allele- and cell-specific CN profiles consistent?** Compares Unique Profile Counts (UPC) and Unique Profile Size Deviation (UPSD) between callers. It identifies if one caller tends to smooth profiles while the other fragments them, or if they agree on the diversity of unique genomic strings.

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

**Are unique haplotype- and cell-specific CN profiles consistent?** Uses UPC and UPSD to evaluate if two callers identify the same unique sets of haplotype-specific genomic profiles.

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

This section evaluates whether callers produce consistent CN phylogenies in the absence of simulated evolutionary ground truth.

### acCNA phylogeny

**Are allele- and cell-specific CN phylogenies consistent?** Compares the Parsimony Scores (PS) of the phylogenies reconstructed from each caller's predictions. Similar scores suggest consistency in the inferred evolutionary complexity.

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

**Are haplotype- and cell-specific CN phylogenies consistent?** Compares Parsimony Scores derived from phased CN profiles to check if the inferred single-cell evolutionary trees have similar structural complexity.

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

This section evaluates the consistency of the accumulated "burden" of complex hcCNA events (focal, medium, and broad) between two callers in the absence of a ground truth.

### Complex hcCNAs

**Are accumulated complex hcCNAs consistent?** Quantifies the overlap of detected complex events using the Jaccard Similarity (JS) index. Within the intersecting genomic regions, it measures the numerical agreement of the copy number strings using RMSE, SCC, and the exact match rate (ACC).

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

By default, `complex_hccnas()` annotates the input CN profiles and compares complex hcCNA calls between tools. Use `genome="hg38"` (default) or `genome="hg19"` to match the reference genome of your data. The default annotation threshold is `0.95`; the example above uses `0.8`.

To save the annotated events for each tool, set `save_annotated_df=True`. The files are written to `{output_dir}/{outprefix}annotated/` as `{tool}_annotated.csv`. You can also use `chromosomes=["chr1", "chr2"]` to limit the analysis to selected chromosomes; both `"chr1"` and `"1"` formats are accepted. By default, all chromosomes are included.

#### Use existing annotated CNA files

If your CN profiles have already been annotated, provide the annotation tables directly to avoid repeating this step:

`exp/` does not contain annotated event CSVs. First generate independent annotations from its `haplotype_combined.csv` files using the preceding raw-input call with `save_annotated_df=True`, or the manual Parser method. Then define `annotated_files` as a mapping from `"Alleloscope"` and `"CHISEL"` to those actual output files. With the preceding call's outprefix, saved files are under `{output_dir}/complex_hccnas_0.8annotated/`; they are generated outputs, not files supplied in `exp/`. Only after these files exist, use this template:

```python
# realbench_runner.complex_hccnas(
#     annotated_files=annotated_files,
#     outprefix="complex_hccnas_reused_",
# )
```

Each tool requires one combined CSV file with the columns `cell`, `Chrom`, `Start`, `End`, `value`, `Segment`, `size`, and `type`. Coordinates must be positive, 1-based, and inclusive. The files should contain the complete independently annotated events, rather than only focal events or events derived from pairwise-aligned data. You can also pass a list of file paths through `annotated_files`, together with matching `tool_names`. Do not use `tool_cna_files` and `annotated_files` in the same call.

The comparison aligns event intervals using shared genomic boundaries while retaining each tool's original event types, sizes, and CN values. Regions without an annotated event are not automatically classified as neutral. The resulting `overlap.csv` and `metrics.csv` report pairwise agreement across shared cell/region entries and event categories; these metrics are not base-pair-weighted measures of event overlap.

For instructions on generating annotated CNA tables, see [Parser annotation](../parsers/index.md).

### Mirrored-subclonal-CNAs

**Are accumulated Mirrored-Subclonal-CNAs consistent?** Annotates mirrored events across the population and evaluates whether both callers agree on these complex evolutionary outcomes by comparing concatenated copy number vectors across the relevant segments and cells.

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

------

## CN Phasing

This section focuses on the consistency of haplotype phasing between callers. It determines if two different algorithms assign the same parental phase to heterozygous copy number states.

### Global CN phasing

**Are CN values consistently phased?** Quantifies the pairwise phasing agreement by calculating the mismatch error between the callers' predicted phase strings. This is conducted in **shared-heterozygous mode**, meaning only loci identified as informative by both callers are compared. A mismatch error near 0 indicates that the two algorithms are in high agreement regarding the global and local phase.

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


## Focal genes and pathways analysis

**Are focal hcCNA-associated genes and pathways consistent across callers and resolutions?** Annotates genes overlapping focal hcCNA regions, performs over-representation analysis (ORA), and calculates gene/pathway Jaccard similarity between callers and resolutions. 

### Prepare gene and pathway annotations

This analysis requires annotated hcCNA events, a gene annotation (GTF) file, and pathway gene sets. HCBench includes three pathway libraries: `MSigDB_Hallmark_2020`, `Reactome_2022`, and `KEGG_2021_Human`.

For gene annotation, HCBench uses [GENCODE v46 (hg38)](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz) or [GENCODE v46lift37 (hg19)](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/gencode.v46lift37.annotation.gtf.gz). If `gtf_path` is not specified, the appropriate GTF file is downloaded automatically when needed and cached in `~/.cache/hcbench/resources/genes/<genome>/`. Alternatively, provide the path to a local GTF file using `gtf_path`. Make sure `genome` matches the genomic coordinates in your CNA files.

By default, all three bundled pathway libraries are available. Use `libraries` to select specific libraries, `gene_sets` to provide custom GMT files, or `resource_dir` to use a custom resource directory. When using the bundled libraries, you can leave `resource_dir` unset.

### Run focal genes and pathways analysis

`exp/` supplies CNA matrices, not annotated hcCNA event tables. Generate the Alleloscope and CHISEL annotations as described above, then supply two `annotated_inputs` records in `tools` order, each with its actual `sample`, `tool`, `bin_size`, and generated annotation `files`. The sample identity and resolution are not inferred from the directory names and are not invented here. For an offline run, also supply an existing GTF as `gtf_path`; `exp/` has no GTF. Omitting that argument retains the automatic GTF download behavior described above. The following template requires these inputs to be prepared first:

```python
# Focal genes and pathways analysis (after preparing the missing inputs)
# realbench_runner.focal_genes_and_pathways(
#     annotated_inputs=annotated_inputs,
#     genome="hg38",  # must match the input coordinates
#     gtf_path=gtf_path,  # existing GTF; omit to use the automatic default
#     workers=2,
#     # outprefix="focal_genes_and_pathways"
# )
```

Each entry in `annotated_inputs` identifies a sample, a caller, a bin size, and its annotated CNA file(s). You can provide either a combined CSV covering multiple chromosomes or a list of chromosome-specific CSV files. Multiple entries allow comparisons across callers and resolutions.

The analysis reports genes overlapping focal hcCNA events, over-representation analysis (ORA) results, gene- and pathway-level Jaccard similarities, and summary statistics. See the output files for results from each sample, caller, and resolution.
