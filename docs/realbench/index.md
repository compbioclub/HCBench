# RealBench

`RealBench` is an integrated module for evaluating copy number alteration (CNA) results in **real-world datasets**. Since real data often lacks a ground truth, this module focuses on cross-tool consistency, evolutionary plausibility, and statistical validation using allele frequencies (VAF) and read depth (RDR).

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

