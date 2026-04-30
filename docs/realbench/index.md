# HCBench-real

`HCBench-real` is an integrated module for evaluating copy number alteration (CNA) results in **real-world datasets**. Since real data often lacks a ground truth, this module focuses on cross-tool consistency, evolutionary plausibility, and statistical validation using allele frequencies (VAF) and read depth (RDR).

------

## Perspective A: Clone Detection

This perspective evaluates the pairwise consistency of clone architectures inferred by two different callers. It determines whether different algorithms agree on the underlying population structure of a real tumor sample.

- **A1: Are global tumor clones consistent?** 

  Evaluates the agreement of clustering labels between two callers using AMI and ARI. These metrics quantify how similar the overall cell groupings are across the entire population.

- **A2: Are rare and dominant clones consistent?** 

  Assesses consistency stratified by predicted clone size. Using Clone Size Deviation (CSD), it examines whether a cell assigned to a cluster of size $n$ by one caller is assigned to a cluster of similar size by the other, identifying if callers disagree on the granularity of specific subclones.



```Python
from hcbench.realbench import RealBench


output_dir = f"/mnt/cbc_adam/public/workspace/HCDSIM/hcbench/input/5Mb/{dataset_name}"
realbench_runner = RealBench(output_dir=f"{output_dir}/new_output")

tools = ["CHISEL", "CNRein", "Alleloscope", "SIGNALS", "SEACON"]
#A1
realbench_runner.clusterConsistency(
        tool_cluster_files=[
            f"{output_dir}/chisel/clusters.csv",
            f"{output_dir}/Alleloscope/clusters.csv",
            f"{output_dir}/signals/clusters.csv"],
        tool_names=["CHISEL", "Alleloscope","SIGNALS"],
)
#A2
realbench_runner.cloneSizebycluster(
        tool_cluster_files=[
            f"{output_dir}/chisel/clusters.csv",
            f"{output_dir}/Alleloscope/clusters.csv",
            f"{output_dir}/signals/clusters.csv"],
        tool_names=["CHISEL", "Alleloscope","SIGNALS"],
)


```

------

## Perspective B: acCNA Detection 

**This perspective assesses the consistency of allele-specific copy number alterations (acCNAs) between callers. In real data scenarios, clone-level comparisons are bypassed in favor of direct cell-level profile matching.**

- **B1: Are unique allele- and cell-specific CN profiles consistent?** Compares Unique Profile Counts (UPC) and Unique Profile Size Deviation (UPSD) between callers. It identifies if one caller tends to smooth profiles while the other fragments them, or if they agree on the diversity of unique genomic strings.
- **B2: Are allele- and cell-specific CN phylogenies consistent?** Compares the Parsimony Scores (PS) of the phylogenies reconstructed from each caller's predictions. Similar scores suggest consistency in the inferred evolutionary complexity.
- **B3: Are allele-, cell-, and locus-specific CN values consistent?** Quantifies the direct numerical agreement of copy number values across all shared cells and loci using RMSE, SCC, and ACC.
- **B4: Are allele-, cell-, and locus-specific CN loss, gain, and neutral consistent?** Uses a binary classification framework to measure concordance. By treating one caller as a reference, metrics like ACC and the Kappa score are used to determine how often the two tools agree on the presence of gains, losses, or neutral states.

```python
# B1
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

# B2
realbench_runner.dolazactree(
        tool_cna_files=[
            f"{output_dir}/chisel/cell_level/minor_major.csv",
            f"{output_dir}/CNRein/minor_major.csv",
            f"{output_dir}/Alleloscope/minor_major.csv",
            f"{output_dir}/signals/minor_major.csv",
            f"{output_dir}/SEACON/minor_major.csv"],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
    )
#B3
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
#B4
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

------

## Perspective C: Single hcCNA Detection 

**This perspective focuses on the consistency of haplotype-specific copy number alterations (hcCNAs). It mirrors the acCNA consistency checks but requires the callers to agree on the specific phasing of the copy number states.**

- **C1: Are unique haplotype- and cell-specific CN profiles consistent?** Uses UPC and UPSD to evaluate if two callers identify the same unique sets of haplotype-specific genomic profiles.
- **C2: Are haplotype- and cell-specific CN phylogenies consistent?** Compares Parsimony Scores derived from phased CN profiles to check if the inferred single-cell evolutionary trees have similar structural complexity.
- **C3: Are haplotype-, cell-, and locus-specific CN values consistent?** Measures numerical consistency for phased copy numbers (e.g., matching 1|2 vs 1|2) using RMSE, SCC, and the exact match rate (ACC).
- **C4: Are haplotype-, cell-, and locus-specific CN loss, gain, and neutral consistent?** Evaluates the reliability of phased CN calling by measuring the concordance (ACC and Kappa score) between callers when categorizing segments into gain, loss, or neutral states.

```Python
# C1
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

# C2
realbench_runner.dolazactree(
        tool_cna_files=[
            f"{output_dir}/chisel/cell_level/haplotype_combined.csv",
            f"{output_dir}/CNRein/haplotype_combined.csv",
            f"{output_dir}/Alleloscope/haplotype_combined.csv",
            f"{output_dir}/signals/haplotype_combined.csv",
            f"{output_dir}/SEACON/haplotype_combined.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
    )

#C3
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
#C4
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

Here are the final three perspectives (D, E, and F) for the real-world data analysis (HCBench-real) in Markdown format.

------

## Perspective D: RDR & VAF Validation 

**This perspective assesses the physical plausibility of the predicted copy numbers by validating them against raw sequencing signals, such as Read Depth Ratio (RDR) and Variant Allele Frequency (VAF).**

- **D1: Are total CNs consistent with observed read depths?** Evaluates the goodness of fit between the inferred total copy number and the observed read-depth signal using the Read Depth L1 Error. Statistically compares callers using paired t-tests; a lower L1 error indicates that the predicted states more accurately reflect the raw data counts.
- **D2: Are major and minor CN distributions consistent with observed VAFs?** Measures how well the predicted major and minor alleles align with truncal SNV frequencies (or BAFs) using log-likelihood under a binomial model. A higher log-likelihood ratio (LLR) indicates a superior fit to the observed molecular evidence.

```python
#D1
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
#D2
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

## Perspective E: Accumulated Complex hcCNA Detection

**This perspective evaluates the consistency of the accumulated "burden" of complex hcCNA events (focal, medium, and broad) between two callers in the absence of a ground truth.**

- **E1: Are accumulated complex hcCNAs consistent?** Quantifies the overlap of detected complex events using the Jaccard Similarity (JS) index. Within the intersecting genomic regions, it measures the numerical agreement of the copy number strings using RMSE, SCC, and the exact match rate (ACC).
- **E2: Are accumulated Mirrored-Subclonal-CNAs consistent?** Annotates mirrored events across the population and evaluates whether both callers agree on these complex evolutionary outcomes by comparing concatenated copy number vectors across the relevant segments and cells.



```python
#E1
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
# E2
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

## Perspective F: Global CN Phasing

**This perspective focuses on the consistency of haplotype phasing between callers. It determines if two different algorithms assign the same parental phase to heterozygous copy number states.**

- **F1: Are CN values consistently phased?** Quantifies the pairwise phasing agreement by calculating the mismatch error between the callers' predicted phase strings. This is conducted in **shared-heterozygous mode**, meaning only loci identified as informative by both callers are compared. A mismatch error near 0 indicates that the two algorithms are in high agreement regarding the global and local phase.



```python
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

















## Overview of Real-Data Submodules

| **Submodule**         | **Description**                                              | **Core Metrics**                         |
| --------------------- | ------------------------------------------------------------ | ---------------------------------------- |
| `clusterConsistency`  | Evaluates the agreement of cluster assignments between different tools. | AMI, ARI                                 |
| `cellprofile`         | Assesses heterogeneity by counting unique genomic profiles identified by each tool. | Unique Profile Count                     |
| `dolazactree`         | Calculates parsimony scores; lower scores indicate more plausible evolutionary paths. | Parsimony Score                          |
| `cndetect`            | Evaluates copy number (CN) detection accuracy                | RMSE, ACC, SCC                           |
| `cnclass`             | Calculates CN state classification metrics                   | AUROC, AUPRC, ACC, Precision, Recall, F1 |
| `segmentation`        | Evaluates the consistency of genomic segmentation boundaries across tools. | Boundary Similarity                      |
| `get_cell_overlap`    | Identifies and analyzes the set of common cells successfully processed by all tools. | Overlap Count / Cell List                |
| `calLLR`              | Uses SNV data to calculate Log-Likelihood Ratios for validating haplotype predictions. | LLR Score                                |
| `findVAFDistribution` | Evaluates VAF distribution consistency for specific copy number states. | VAF Density                              |
| `rddetect`            | Checks if predicted copy numbers align with observed Read Depth Ratios. | RDR Consistency                          |
| `hcPhasing`           |                                                              |                                          |

------

## Example Usage

The following Python script demonstrates a full benchmark run on a real-world dataset, including the newly added segmentation and overlap modules.

Python

```
from hcbench.realbench import RealBench


output_dir = f"/mnt/cbc_adam/public/workspace/HCDSIM/hcbench/input/5Mb/{dataset_name}"
realbench_runner = RealBench(output_dir=f"{output_dir}/new_output")

# Tool definitions
tools = ["CHISEL", "CNRein", "Alleloscope", "SIGNALS", "SEACON"]

realbench_runner.clusterConsistency(
        tool_cluster_files=[
            f"{output_dir}/chisel/clusters.csv",
            f"{output_dir}/Alleloscope/clusters.csv",
            f"{output_dir}/signals/clusters.csv"],
        tool_names=["CHISEL", "Alleloscope","SIGNALS"],
)
    
realbench_runner.cloneSizebycluster(
        tool_cluster_files=[
            f"{output_dir}/chisel/clusters.csv",
            f"{output_dir}/Alleloscope/clusters.csv",
            f"{output_dir}/signals/clusters.csv"],
        tool_names=["CHISEL", "Alleloscope","SIGNALS"],
)


# Count unique cell profiles (Haplotype-combined)
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

# Compute Parsimony Scores
realbench_runner.dolazactree(
        tool_cna_files=[
            f"{output_dir}/chisel/cell_level/haplotype_combined.csv",
            f"{output_dir}/CNRein/haplotype_combined.csv",
            f"{output_dir}/Alleloscope/haplotype_combined.csv",
            f"{output_dir}/signals/haplotype_combined.csv",
            f"{output_dir}/SEACON/haplotype_combined.csv",],
        tool_names=["CHISEL", "CNRein","Alleloscope", "SIGNALS","SEACON"],
    )
    
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

# ============================================================
# 3) Statistical Validation (VAF & RDR)
# ============================================================

# Calculate Log-Likelihood Ratio (LLR) for predicted haplotypes
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

# Check RDR (Read Depth Ratio) detection accuracy
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

realbench_runner.get_cell_overlap(
    tool_cna_files=[
        f"{output_dir}/chisel/cell_level/haplotype_combined.csv",
        f"{output_dir}/CNRein/haplotype_combined.csv",
        f"{output_dir}/Alleloscope/haplotype_combined.csv",
        f"{output_dir}/signals/haplotype_combined.csv",
        f"{output_dir}/SEACON/haplotype_combined.csv",
    ],
    tool_names=tools,
)

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
    mode = "homozygous-inclusive-all",
    # outprefix = "hcPhasing_homorozygous-included"
)

```



------

