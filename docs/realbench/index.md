# HCBench-real

`HCBench-real` is an integrated module for evaluating copy number alteration (CNA) results in **real-world datasets**. Since real data often lacks a ground truth, this module focuses on pairwise cross-tool consistency and agreement with observed sequencing signals. The evaluation is organized around five evaluation questions: CNA detection, clone detection, evolution tracking, accumulated complex hcCNA detection, and CN phasing.

## Setup

```python
from hcbench.realbench import RealBench

output_dir = f"/mnt/cbc_adam/public/workspace/HCDSIM/hcbench/input/5Mb/{dataset_name}"
realbench_runner = RealBench(output_dir=f"{output_dir}/new_output")

tools = ["CHISEL", "CNRein", "Alleloscope", "SIGNALS", "SEACON"]
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
            f"{output_dir}/chisel/cell_level/minor_major.csv",
            f"{output_dir}/CNRein/minor_major.csv",
            f"{output_dir}/Alleloscope/minor_major.csv",
            f"{output_dir}/signals/minor_major.csv",
            f"{output_dir}/SEACON/minor_major.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
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
    tool_hap1_cna_files = [
        f"{output_dir}/chisel/cell_level/minor.csv",
            f"{output_dir}/Alleloscope/minor.csv",
            f"{output_dir}/CNRein/minor.csv",
            f"{output_dir}/SEACON/minor.csv",
            f"{output_dir}/signals/minor.csv"],
    tool_hap2_cna_files = [
        f"{output_dir}/chisel/cell_level/major.csv",
            f"{output_dir}/Alleloscope/major.csv",
            f"{output_dir}/CNRein/major.csv",
            f"{output_dir}/SEACON/major.csv",
            f"{output_dir}/signals/major.csv"],
    tool_names = ["CHISEL", "Alleloscope","CNRein", "SEACON","SIGNALS"], 
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
            f"{output_dir}/chisel/cell_level/haplotype_combined.csv",
            f"{output_dir}/CNRein/haplotype_combined.csv",
            f"{output_dir}/Alleloscope/haplotype_combined.csv",
            f"{output_dir}/signals/haplotype_combined.csv",
            f"{output_dir}/SEACON/haplotype_combined.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
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
    tool_hap1_cna_files = [
        f"{output_dir}/chisel/cell_level/haplotype_1.csv",
        f"{output_dir}/Alleloscope/haplotype_1.csv",
        f"{output_dir}/CNRein/haplotype_1.csv",
        f"{output_dir}/SEACON/haplotype_1.csv",
        f"{output_dir}/signals/haplotype_1.csv"],
    tool_hap2_cna_files = [
        f"{output_dir}/chisel/cell_level/haplotype_2.csv",
        f"{output_dir}/Alleloscope/haplotype_2.csv",
        f"{output_dir}/CNRein/haplotype_2.csv",
        f"{output_dir}/SEACON/haplotype_2.csv",
        f"{output_dir}/signals/haplotype_2.csv"],
    tool_names = ["CHISEL", "Alleloscope","CNRein", "SEACON","SIGNALS"], 
    outfile="cna_states_hccna.csv",
)
```

### Total CN vs. RDR

**Are total CNs consistent with observed read depths?** Evaluates the goodness of fit between the inferred total copy number and the observed read-depth signal using the Read Depth L1 Error. Statistically compares callers using paired t-tests; a lower L1 error indicates that the predicted states more accurately reflect the raw data counts.

```python
# Total CN vs. RDR
realbench_runner.total_cn_vs_rdr(
        bin_count_files=[
            f"{output_dir}/chisel/bin_counts.csv",
            f"{output_dir}/Alleloscope/bin_counts.csv",
            f"{output_dir}/SEACON/bin_counts.csv",
            f"{output_dir}/CNRein/bin_rdr.csv",
            f"{output_dir}/signals/bin_counts.csv",
            ],
        tool_cna_files=[
            f"{output_dir}/chisel/cell_level/haplotype_combined.csv",
            f"{output_dir}/Alleloscope/haplotype_combined.csv",
            f"{output_dir}/SEACON/haplotype_combined.csv",
            f"{output_dir}/CNRein/haplotype_combined.csv",
            f"{output_dir}/signals/haplotype_combined.csv"],
        tool_names=["CHISEL","Alleloscope", "SEACON","CNRein",'SIGNALS'],
)
```

### Major/minor CN vs. VAF

**Are major and minor CN distributions consistent with observed VAFs?** Measures how well the predicted major and minor alleles align with truncal SNV frequencies (or BAFs) using log-likelihood under a binomial model. A higher log-likelihood ratio (LLR) indicates a superior fit to the observed molecular evidence.

```python
# Major/minor CN vs. VAF
realbench_runner.major_minor_cn_vs_vaf(
    tool_hap1_files = [
        f"{output_dir}/chisel/cell_level/haplotype_1.csv",
        f"{output_dir}/Alleloscope/haplotype_1.csv",
        f"{output_dir}/SEACON/haplotype_1.csv",
        f"{output_dir}/CNRein/haplotype_1.csv",
        f"{output_dir}/signals/haplotype_1.csv"],
    tool_hap2_files = [
        f"{output_dir}/chisel/cell_level/haplotype_2.csv",
        f"{output_dir}/Alleloscope/haplotype_2.csv",
        f"{output_dir}/SEACON/haplotype_2.csv",
        f"{output_dir}/CNRein/haplotype_2.csv",
        f"{output_dir}/signals/haplotype_2.csv"],
    tool_names = ["CHISEL","Alleloscope", "SEACON","CNRein",'SIGNALS'],
    snv_paths  =[
        f"{output_dir}/chisel/VAF/",
        f"{output_dir}/Alleloscope/VAF/",
        f"{output_dir}/SEACON/VAF/",
        f"{output_dir}/CNRein/VAF/",
        f"{output_dir}/signals/VAF/"],
    cell_lists  = cell_lists,
    variant_pos_dfs = variant_pos_dfs,
)
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
            f"{output_dir}/chisel/clusters.csv",
            f"{output_dir}/Alleloscope/clusters.csv",
            f"{output_dir}/signals/clusters.csv"],
        tool_names=["CHISEL", "Alleloscope","SIGNALS"],
)
```

### Rare and dominant clones

**Are rare and dominant clones consistent?**

Assesses consistency stratified by predicted clone size. Using Clone Size Deviation (CSD), it examines whether a cell assigned to a cluster of size $n$ by one caller is assigned to a cluster of similar size by the other, identifying if callers disagree on the granularity of specific subclones.

```python
# Rare and dominant clones
realbench_runner.rare_and_dominant_clones(
        tool_cluster_files=[
            f"{output_dir}/chisel/clusters.csv",
            f"{output_dir}/Alleloscope/clusters.csv",
            f"{output_dir}/signals/clusters.csv"],
        tool_names=["CHISEL", "Alleloscope","SIGNALS"],
)
```

### Unique acCNA profiles

**Are unique allele- and cell-specific CN profiles consistent?** Compares Unique Profile Counts (UPC) and Unique Profile Size Deviation (UPSD) between callers. It identifies if one caller tends to smooth profiles while the other fragments them, or if they agree on the diversity of unique genomic strings.

```python
# Unique acCNA profiles
realbench_runner.unique_cna_profile_counts(
        tool_cna_files=[
             f"{output_dir}/chisel/cell_level/minor_major.csv",
            f"{output_dir}/CNRein/minor_major.csv",
            f"{output_dir}/Alleloscope/minor_major.csv",
            f"{output_dir}/signals/minor_major.csv",
            f"{output_dir}/SEACON/minor_major.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
    outfile="unique_cna_profile_counts_accna.csv",
    )

realbench_runner.unique_cna_profile_sizes(
        tool_cna_files=[
             f"{output_dir}/chisel/cell_level/minor_major.csv",
            f"{output_dir}/CNRein/minor_major.csv",
            f"{output_dir}/Alleloscope/minor_major.csv",
            f"{output_dir}/signals/minor_major.csv",
            f"{output_dir}/SEACON/minor_major.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
    outfile="unique_cna_profile_sizes_accna.csv",
    )
```

### Unique hcCNA profiles

**Are unique haplotype- and cell-specific CN profiles consistent?** Uses UPC and UPSD to evaluate if two callers identify the same unique sets of haplotype-specific genomic profiles.

```python
# Unique hcCNA profiles
realbench_runner.unique_cna_profile_counts(
        tool_cna_files=[
            f"{output_dir}/chisel/cell_level/haplotype_combined.csv",
            f"{output_dir}/CNRein/haplotype_combined.csv",
            f"{output_dir}/Alleloscope/haplotype_combined.csv",
            f"{output_dir}/signals/haplotype_combined.csv",
            f"{output_dir}/SEACON/haplotype_combined.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
    outfile="unique_cna_profile_counts_hccna.csv",
    )

realbench_runner.unique_cna_profile_sizes(
        tool_cna_files=[
            f"{output_dir}/chisel/cell_level/haplotype_combined.csv",
            f"{output_dir}/CNRein/haplotype_combined.csv",
            f"{output_dir}/Alleloscope/haplotype_combined.csv",
            f"{output_dir}/signals/haplotype_combined.csv",
            f"{output_dir}/SEACON/haplotype_combined.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
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
            f"{output_dir}/chisel/cell_level/minor_major.csv",
            f"{output_dir}/CNRein/minor_major.csv",
            f"{output_dir}/Alleloscope/minor_major.csv",
            f"{output_dir}/signals/minor_major.csv",
            f"{output_dir}/SEACON/minor_major.csv"],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
    outfile="cna_phylogeny_accna.csv",
    )
```

### hcCNA phylogeny

**Are haplotype- and cell-specific CN phylogenies consistent?** Compares Parsimony Scores derived from phased CN profiles to check if the inferred single-cell evolutionary trees have similar structural complexity.

```python
# hcCNA phylogeny
realbench_runner.cna_phylogeny(
        tool_cna_files=[
            f"{output_dir}/chisel/cell_level/haplotype_combined.csv",
            f"{output_dir}/CNRein/haplotype_combined.csv",
            f"{output_dir}/Alleloscope/haplotype_combined.csv",
            f"{output_dir}/signals/haplotype_combined.csv",
            f"{output_dir}/SEACON/haplotype_combined.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
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
        f"{output_dir}/chisel/cell_level/haplotype_combined.csv",
        f"{output_dir}/CNRein/haplotype_combined.csv",
        f"{output_dir}/Alleloscope/haplotype_combined.csv",
        f"{output_dir}/signals/haplotype_combined.csv",
        f"{output_dir}/SEACON/haplotype_combined.csv",
    ],
    tool_names=tools,
    threshold=0.8,
    outprefix="complex_hccnas_0.8",
)
```

`complex_hccnas` also accepts `save_annotated_df=False` and
`chromosomes=None`. `None` uses all input chromosomes; for example,
`chromosomes=["chr1", "chr2"]` restricts annotation and evaluation to those
chromosomes (`"1"` and `"2"` are also accepted). An empty list selects none.

With `save_annotated_df=True`, each tool's annotation is saved under
`{output_dir}/{outprefix}annotated/`. The combined file is named
`{tool}_annotated.csv`. If a chromosome list is supplied, additional files
`{tool}_annotated_chr1.csv`, etc. are written for each requested chromosome.
Chromosomes with no annotated events produce a header-only file. Per-chromosome
files are subsets of the combined annotation, not independent re-annotations.
WGD classification retains the existing whole-genome reference denominator.

For real-data comparisons, files are placed in a `{tool1}__{tool2}/`
subdirectory because annotation uses pairwise-aligned bins. Both tools' tables
are saved for each comparison; they are not overwritten by subsequent pairs.
With fewer than two tools, there are no pairwise annotation tables to save.

### Mirrored-subclonal-CNAs

**Are accumulated Mirrored-Subclonal-CNAs consistent?** Annotates mirrored events across the population and evaluates whether both callers agree on these complex evolutionary outcomes by comparing concatenated copy number vectors across the relevant segments and cells.

```python
# Mirrored-subclonal-CNAs
realbench_runner.mirrored_subclonal_cnas(
        tool_cna_files=[
            f"{output_dir}/chisel/cell_level/haplotype_combined.csv",
            f"{output_dir}/CNRein/haplotype_combined.csv",
            f"{output_dir}/Alleloscope/haplotype_combined.csv",
            f"{output_dir}/signals/haplotype_combined.csv",
            f"{output_dir}/SEACON/haplotype_combined.csv",],
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
    tool_hap1_cna_files = [
        f"{output_dir}/chisel/cell_level/haplotype_1.csv",
        f"{output_dir}/Alleloscope/haplotype_1.csv",
        f"{output_dir}/SEACON/haplotype_1.csv",
        f"{output_dir}/CNRein/haplotype_1.csv",
        f"{output_dir}/signals/haplotype_1.csv"],
    tool_hap2_cna_files = [
        f"{output_dir}/chisel/cell_level/haplotype_2.csv",
        f"{output_dir}/Alleloscope/haplotype_2.csv",
        f"{output_dir}/SEACON/haplotype_2.csv",
        f"{output_dir}/CNRein/haplotype_2.csv",
        f"{output_dir}/signals/haplotype_2.csv"],
    tool_names = ["CHISEL","Alleloscope", "SEACON","CNRein",'SIGNALS'],
    mode = "heterozygous-only",
    # outprefix = "global_cn_phasing_heterozygous-only"
)
```


## Focal genes and pathways analysis

**Are focal hcCNA-associated genes and pathways consistent across callers and resolutions?** Annotates genes overlapping focal hcCNA regions, performs over-representation analysis (ORA), and calculates gene/pathway Jaccard similarity between callers and resolutions. 

### Prepare the GTF resource

The three GMT pathway libraries (`MSigDB_Hallmark_2020`, `Reactome_2022`, `KEGG_2021_Human`) are bundled in `hcbench/lib/pathways/` and included in the pip package. GTF is not bundled. If `gtf_path` is omitted, HCBench uses a cached default or automatically downloads [GENCODE v46 for hg38](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz) or [v46lift37 for hg19](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/gencode.v46lift37.annotation.gtf.gz). An explicit `gtf_path` uses the existing file without downloading; no release parameter is required. The CNA coordinates must match the selected genome.

The example below uses an explicit GTF path; bundled GMTs are selected automatically. You can still select a subset with `libraries`, supply exact custom GMT paths with `gene_sets`, or use `resource_dir` for a custom GTF/GMT directory. An explicit `resource_dir` selects its `pathways/` directory instead of the bundled GMTs; leave it unset when using the bundled libraries with an explicit GTF path.

### Run focal genes and pathways analysis

```python
# Focal genes and pathways analysis
realbench_runner.focal_genes_and_pathways(
    annotated_inputs=[
        {
            "sample": "sample1",
            "tool": "CHISEL",
            "bin_size": "100kb",
            "files": f"{output_dir}/CHISEL_annotated.csv"
        },
        {
            "sample": "sample1",
            "tool": "SEACON",
            "bin_size": "100kb",
            "files": [
                f"{output_dir}/SEACON_chr1.csv",
                f"{output_dir}/SEACON_chr2.csv"
            ]
        }
    ],
    genome="hg38",
    gtf_path="~/hcbench_resources/genes/hg38/gencode.v46.annotation.gtf.gz",
    workers=2,
    # outprefix="focal_genes_and_pathways"
)
```

With no `gtf_path`, the GTF cache defaults to `~/.cache/hcbench/resources/genes/<genome>/`; the version and path are printed when reused, and the completed path is printed after download. Download failures remove the incomplete temporary file and raise an error. The hg19 default is now v46lift37 rather than v19, so changing annotation resources can change results.

The results include focal gene annotations, ORA enrichment, gene/pathway Jaccard comparisons, and summary statistics. Both combined multi-chromosome CSV files and separate chromosome files are supported.

See the [detailed documentation](focal_genes_and_pathways.md) for additional parameters and output descriptions.