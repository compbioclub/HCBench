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
realbench_runner.cndetect(
        tool_cna_files=[
            f"{output_dir}/chisel/cell_level/minor_major.csv",
            f"{output_dir}/CNRein/minor_major.csv",
            f"{output_dir}/Alleloscope/minor_major.csv",
            f"{output_dir}/signals/minor_major.csv",
            f"{output_dir}/SEACON/minor_major.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
        haplotype = "combined",
        outfile_prefix= "bin_level"
    )
```

### acCNA states

**Are allele-, cell-, and locus-specific CN loss, gain, and neutral states consistent?**

Uses a binary classification framework to measure concordance. By treating one caller as a reference, metrics like ACC and the Kappa score are used to determine how often the two tools agree on the presence of gains, losses, or neutral states.

```python
# acCNA states
realbench_runner.cnclass(
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
)
```

### hcCNA values

**Are haplotype-, cell-, and locus-specific CN values consistent?**

Measures numerical consistency for phased copy numbers (e.g., matching 1|2 vs 1|2) using RMSE, SCC, and the exact match rate (ACC).

```python
# hcCNA values
realbench_runner.cndetect(
        tool_cna_files=[
            f"{output_dir}/chisel/cell_level/haplotype_combined.csv",
            f"{output_dir}/CNRein/haplotype_combined.csv",
            f"{output_dir}/Alleloscope/haplotype_combined.csv",
            f"{output_dir}/signals/haplotype_combined.csv",
            f"{output_dir}/SEACON/haplotype_combined.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
        haplotype = "combined",
        outfile_prefix= "bin_level"
    )
```

### hcCNA states

**Are haplotype-, cell-, and locus-specific CN loss, gain, and neutral states consistent?**

Evaluates the reliability of phased CN calling by measuring the concordance (ACC and Kappa score) between callers when categorizing segments into gain, loss, or neutral states.

```python
# hcCNA states
realbench_runner.cnclass(
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
)
```

### Total CN vs. RDR

**Are total CNs consistent with observed read depths?** Evaluates the goodness of fit between the inferred total copy number and the observed read-depth signal using the Read Depth L1 Error. Statistically compares callers using paired t-tests; a lower L1 error indicates that the predicted states more accurately reflect the raw data counts.

```python
# Total CN vs. RDR
realbench_runner.rddetect(
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
realbench_runner.calLLR(
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
realbench_runner.clusterConsistency(
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
realbench_runner.cloneSizebycluster(
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
realbench_runner.cellprofile(
        tool_cna_files=[
             f"{output_dir}/chisel/cell_level/minor_major.csv",
            f"{output_dir}/CNRein/minor_major.csv",
            f"{output_dir}/Alleloscope/minor_major.csv",
            f"{output_dir}/signals/minor_major.csv",
            f"{output_dir}/SEACON/minor_major.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
    )

realbench_runner.cloneSizebycellprofile(
        tool_cna_files=[
             f"{output_dir}/chisel/cell_level/minor_major.csv",
            f"{output_dir}/CNRein/minor_major.csv",
            f"{output_dir}/Alleloscope/minor_major.csv",
            f"{output_dir}/signals/minor_major.csv",
            f"{output_dir}/SEACON/minor_major.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
    )
```

### Unique hcCNA profiles

**Are unique haplotype- and cell-specific CN profiles consistent?** Uses UPC and UPSD to evaluate if two callers identify the same unique sets of haplotype-specific genomic profiles.

```python
# Unique hcCNA profiles
realbench_runner.cellprofile(
        tool_cna_files=[
            f"{output_dir}/chisel/cell_level/haplotype_combined.csv",
            f"{output_dir}/CNRein/haplotype_combined.csv",
            f"{output_dir}/Alleloscope/haplotype_combined.csv",
            f"{output_dir}/signals/haplotype_combined.csv",
            f"{output_dir}/SEACON/haplotype_combined.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
    )

realbench_runner.cloneSizebycellprofile(
        tool_cna_files=[
            f"{output_dir}/chisel/cell_level/haplotype_combined.csv",
            f"{output_dir}/CNRein/haplotype_combined.csv",
            f"{output_dir}/Alleloscope/haplotype_combined.csv",
            f"{output_dir}/signals/haplotype_combined.csv",
            f"{output_dir}/SEACON/haplotype_combined.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
    )
```

------

## Evolution Tracking

This section evaluates whether callers produce consistent CN phylogenies in the absence of simulated evolutionary ground truth.

### acCNA phylogeny

**Are allele- and cell-specific CN phylogenies consistent?** Compares the Parsimony Scores (PS) of the phylogenies reconstructed from each caller's predictions. Similar scores suggest consistency in the inferred evolutionary complexity.

```python
# acCNA phylogeny
realbench_runner.dolazactree(
        tool_cna_files=[
            f"{output_dir}/chisel/cell_level/minor_major.csv",
            f"{output_dir}/CNRein/minor_major.csv",
            f"{output_dir}/Alleloscope/minor_major.csv",
            f"{output_dir}/signals/minor_major.csv",
            f"{output_dir}/SEACON/minor_major.csv"],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
    )
```

### hcCNA phylogeny

**Are haplotype- and cell-specific CN phylogenies consistent?** Compares Parsimony Scores derived from phased CN profiles to check if the inferred single-cell evolutionary trees have similar structural complexity.

```python
# hcCNA phylogeny
realbench_runner.dolazactree(
        tool_cna_files=[
            f"{output_dir}/chisel/cell_level/haplotype_combined.csv",
            f"{output_dir}/CNRein/haplotype_combined.csv",
            f"{output_dir}/Alleloscope/haplotype_combined.csv",
            f"{output_dir}/signals/haplotype_combined.csv",
            f"{output_dir}/SEACON/haplotype_combined.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
    )
```

------

## Accumulated Complex hcCNA Detection

This section evaluates the consistency of the accumulated "burden" of complex hcCNA events (focal, medium, and broad) between two callers in the absence of a ground truth.

### Complex hcCNAs

**Are accumulated complex hcCNAs consistent?** Quantifies the overlap of detected complex events using the Jaccard Similarity (JS) index. Within the intersecting genomic regions, it measures the numerical agreement of the copy number strings using RMSE, SCC, and the exact match rate (ACC).

```python
# Complex hcCNAs
realbench_runner.segmentation(
    tool_cna_files=[
        f"{output_dir}/chisel/cell_level/haplotype_combined.csv",
        f"{output_dir}/CNRein/haplotype_combined.csv",
        f"{output_dir}/Alleloscope/haplotype_combined.csv",
        f"{output_dir}/signals/haplotype_combined.csv",
        f"{output_dir}/SEACON/haplotype_combined.csv",
    ],
    tool_names=tools,
    threshold=0.8,
    outprefix="segmentation_0.8",
)
```

### Mirrored-subclonal-CNAs

**Are accumulated Mirrored-Subclonal-CNAs consistent?** Annotates mirrored events across the population and evaluates whether both callers agree on these complex evolutionary outcomes by comparing concatenated copy number vectors across the relevant segments and cells.

```python
# Mirrored-subclonal-CNAs
 realbench_runner.get_mirrored_subclonal(
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
realbench_runner.hcPhasing(
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
    # outprefix = "hcPhasing_heterozygous-only"
)
```

