# HCBench-sim

`HCBench-sim` evaluate inferred hcCNA profiles against the hcCNA-bench GT. The evaluation comprised three baseline perspectives and three advanced perspectives.

------

## Perspective A: Clone Detection

This perspective evaluates the accuracy of inferring clonal architecture within a tumor population. It focuses on whether cells are correctly grouped into their ancestral lineages, which is a prerequisite for interpreting copy number alterations.

- **A1: Are global tumor clones correct?**

  Uses global metrics like Adjusted Mutual Information (AMI) and Adjusted Rand Index (ARI) to compare predicted labels against the ground truth. This provides a high-level assessment of how well the overall population structure is captured.

- **A2: Are rare and dominant clones correct?** 

  Evaluates clustering accuracy stratified by clone size using Clone Size Deviation (CSD). It measures whether specific clones (regardless of size) are being merged into larger clusters (positive deviation) or "shattered" into many small groups (negative deviation).

```Python
from hcbench.gtbench.gtbench import GTBench

# All benchmark outputs will be written under this directory
bench = GTBench(output_dir="/path/to/results/gt_output/cell_level")

# A1
bench.clusterConsistency(
    tool_clone_files=[
        "/path/to/signals/clusters.csv",
        "/path/to/seacon/clusters.csv",
    ],
    tool_names=[
        "signals",
        "seacon",
    ],
)

# A2
bench.cloneSizebycluster(
    tool_cluster_files=[
        "/path/to/signals/clusters.csv",
        "/path/to/seacon/clusters.csv",
    ],
    gt_cluster_file="/path/to/gt/clusters.csv",
    tool_names=[
        "signals",
        "seacon",
    ],
    outfile="clone_size_by_cluster.csv",
)
```

------

## Perspective B: acCNA Detection

This perspective evaluates the detection of allele-specific copy number alterations (acCNAs). It addresses the accuracy of unique genomic profiles, phylogenetic reconstruction, and specific copy number values at both the single-cell and clone levels.

- **B1: Are unique allele- and cell-specific CN profiles correct?** 

  Analyzes "Unique Profile Groups" (UPGs) using Unique Profile Count Error (UPCE) and Unique Profile Size Deviation (UPSD). These metrics detect if distinct biological profiles are being over-smoothed (merged) or over-fragmented due to noise.

- **B2: Are allele- and cell-specific CN phylogenies correct?** 

  Applies a maximum parsimony criterion to reconstructed cell lineages. It measures Parsimony Score Deviation (PSD) to see if the inferred evolutionary history is artificially complex (indicating noise) or overly simplified.

- **B3: Are allele-, cell/clone-, and locus-specific CN values correct?** 

  Assesses quantitative accuracy using RMSE and SCC for major and minor alleles. It also uses ACC to measure the percentage of loci where the full allele-specific string is exactly correct.

- **B4: Are allele-, cell/clone-, and locus-specific CN loss, gain, and neutral correct?**

  Treats copy number detection as a binary classification problem (Gain, Neutral, Loss). Performance is measured via a suite of statistical metrics (AUROC, Sensitivity, F1-Score, etc.) to determine the reliability of calling specific genomic states.

```Python
# B1
bench.cellprofile(
    tool_cna_files=[
        "/path/to/gt/cell/minor_major.csv",
        "/path/to/signals/minor_major.csv",
        "/path/to/seacon/minor_major.csv",
    ],
    tool_names=[
        "GT",
        "signals",
        "seacon",
    ],
    outfile="unique_cell_profile_accna.csv",
)

bench.cloneSizebycellprofile(
    tool_cna_files=[
        "/path/to/signals/minor_major.csv",
        "/path/to/seacon/minor_major.csv",
    ],
    gt_cna_file="/path/to/gt/cell/minor_major.csv",
    tool_names=[
        "signals",
        "seacon",
    ],
    outfile="clone_size_by_cell_profile_accna.csv",
)
# B2
bench.dolazactree(
    tool_cna_files=[
        "/path/to/gt/cell/minor_major.csv",
        "/path/to/signals/minor_major.csv",
        "/path/to/seacon/minor_major.csv",
    ],
    tool_names=[
        "GT",
        "signals",
        "seacon",
    ],
    outfile="parsimony_score_accna.csv",
)

#B3
bench.cndetect(
    tool_cna_files=[
        "/path/to/signals/minor_major.csv",
        "/path/to/seacon/minor_major.csv",
    ],
    cna_profile_file="/path/to/gt/cell/minor_major.csv",
    tool_names=[
        "signals",
        "seacon",
    ],
    haplotype="combined",
    outfile="bin_level_results_accna.csv",
)
#B4
bench.cnclass(
    tool_hap1_cna_files=[
        "/path/to/signals/minor.csv",
        "/path/to/seacon/minor.csv",
    ],
    tool_hap2_cna_files=[
        "/path/to/signals/major.csv",
        "/path/to/seacon/major.csv",
    ],
    tool_names=[
        "signals",
        "seacon",
    ],
    profile_hap1_cna_file="/path/to/gt/cell/minor.csv",
    profile_hap2_cna_file="/path/to/gt/cell/major.csv",
    outfile="cnclass_results_accna.csv",
)

```

------



## Perspective C: Single hcCNA Detection

This perspective focuses on haplotype-specific copy number alterations (hcCNAs). It mirrors the logic of acCNA detection but requires higher precision by distinguishing between the specific phased haplotypes (e.g., 3|1 vs. 1|3).

- **C1: Are unique haplotype- and cell-specific CN profiles correct?**

   Quantifies performance based on cells sharing identical haplotype-specific CN strings across the genome using UPCE and UPSD metrics.

- **C2: Are haplotype- and cell-specific CN phylogenies correct?** 

  Calculates PSD using phased CN profiles as input to build the phylogeny, ensuring that the evolutionary relationships account for the specific parent-of-origin alleles.

- **C3: Are haplotype-, cell/clone-, and locus-specific CN values correct?** 

  Calculates the average RMSE and SCC across both individual haplotypes and evaluates the rate of exact matches for the phased CN strings (ACC).

- **C4: Are haplotype-, cell/clone-, and locus-specific CN loss, gain, and neutral correct?** 

  Assesses phased CN calling performance via binary classification across multiple states (Gain/Loss/Neutral), providing a comprehensive statistical breakdown of detection sensitivity and specificity.



```python
# C1
bench.cellprofile(
    tool_cna_files=[
        "/path/to/gt/cell/haplotype_combined.csv",
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    tool_names=[
        "GT",
        "signals",
        "seacon",
    ],
    outfile="unique_cell_profile.csv",  # Output CSV name under output_dir
)

bench.cloneSizebycellprofile(
    tool_cna_files=[
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    gt_cna_file="/path/to/gt/cell/haplotype_combined.csv",
    tool_names=[
        "signals",
        "seacon",
    ],
    outfile="clone_size_by_cell_profile.csv",
)
#C2
bench.dolazactree(
    tool_cna_files=[
        "/path/to/gt/cell/haplotype_combined.csv",
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    tool_names=[
        "GT",
        "signals",
        "seacon",
    ],
    outfile="parsimony_score.csv",
)
#C3
bench.cndetect(
    tool_cna_files=[
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    cna_profile_file="/path/to/gt/cell/haplotype_combined.csv",
    tool_names=[
        "signals",
        "seacon",
    ],
    haplotype="combined",
    outfile="bin_level_results.csv",
)
#C4
bench.cnclass(
    tool_hap1_cna_files=[
        "/path/to/signals/haplotype_1.csv",
        "/path/to/seacon/haplotype_1.csv",
    ],
    tool_hap2_cna_files=[
        "/path/to/signals/haplotype_2.csv",
        "/path/to/seacon/haplotype_2.csv",
    ],
    tool_names=[
        "signals",
        "seacon",
    ],
    profile_hap1_cna_file="/path/to/gt/cell/haplotype_1.csv",
    profile_hap2_cna_file="/path/to/gt/cell/haplotype_2.csv",
    type="hcCNA",                       # Keep only if your cnclass() signature supports it
    outfile="cnclass_results.csv",
)

```



------

## Perspective D: Evolving Complex hcCNA Detection

This perspective evaluates the model's ability to identify the temporal onset and evolutionary tracking of complex haplotype-specific copy number alterations (hcCNAs) across various genomic scales (focal, medium, and broad).

- **D1: Are parent CNs of complex hcCNAs correct at onset?** Assesses whether the model correctly identifies the pre-event (parental) copy number state at the exact moment an alteration arises during evolution, using metrics like ACC and RMSE.
- **D2: Are CN changes of complex hcCNAs correct at onset?** Evaluates the model's accuracy in capturing the magnitude and direction of the copy number transition (e.g., a gain from 1 to 3) from parent to child.
- **D3: Are complex hcCNA correct at onset?** A holistic measure requiring the parent state, child state, and the resulting change to all match the ground truth. This reflects the end-to-end reliability in pinpointing true evolutionary transitions.
- **D4: Are complex hcCNA stable throughout evolution?** Evaluates whether a detected alteration is correctly tracked as a stable, heritable state across all descendant clones following its initial onset.

```Python
# ============================================================
#  Evolution-onset benchmarks: run for focal / medium / broad changes
#    - hccnchange: D2
#    - hccnstable: D4
#    - hconsetacc: D3
#    - hconsetcn: D1
# ============================================================
for size in ["focal", "medium", "broad"]:

    bench.hccnchange(
        tool_hap1_cna_files=[
            "/path/to/signals/haplotype_1.csv",
            "/path/to/seacon/haplotype_1.csv",
        ],
        tool_hap2_cna_files=[
            "/path/to/signals/haplotype_2.csv",
            "/path/to/seacon/haplotype_2.csv",
        ],
        tool_names=[
            "signals",
            "seacon",
        ],
        changes_file=f"/path/to/gt/cell/cell_changes_split_annotated_{size}.csv",
        outfile=f"evolution_onset_CN_Change_{size}.csv",
    )

    bench.hccnstable(
        tool_hap1_cna_files=[
            "/path/to/signals/haplotype_1.csv",
            "/path/to/seacon/haplotype_1.csv",
        ],
        tool_hap2_cna_files=[
            "/path/to/signals/haplotype_2.csv",
            "/path/to/seacon/haplotype_2.csv",
        ],
        tool_names=[
            "signals",
            "seacon",
        ],
        changes_file=f"/path/to/gt/cell/cell_changes_split_annotated_{size}.csv",
        tree_file="/path/to/tree/tree_new.newick",
        outfile=f"evolution_cn_stability_acc_{size}.csv",
    )

    bench.hconsetacc(
        tool_hap1_cna_files=[
            "/path/to/signals/haplotype_1.csv",
            "/path/to/seacon/haplotype_1.csv",
        ],
        tool_hap2_cna_files=[
            "/path/to/signals/haplotype_2.csv",
            "/path/to/seacon/haplotype_2.csv",
        ],
        tool_names=[
            "signals",
            "seacon",
        ],
        changes_file=f"/path/to/gt/cell/cell_changes_split_annotated_{size}.csv",
        outfile=f"evolution_onset_acc_{size}.csv",
    )

    bench.hconsetcn(
        tool_hap1_cna_files=[
            "/path/to/signals/haplotype_1.csv",
            "/path/to/seacon/haplotype_1.csv",
        ],
        tool_hap2_cna_files=[
            "/path/to/signals/haplotype_2.csv",
            "/path/to/seacon/haplotype_2.csv",
        ],
        tool_names=[
            "signals",
            "seacon",
        ],
        changes_file=f"/path/to/gt/cell/cell_changes_split_annotated_{size}.csv",
        outfile=f"evolution_onset_parent_CN_{size}.csv",
    )
```

------

## Perspective E: Accumulated Complex hcCNA Detection

Perspective E assesses the final "burden" of complex hcCNAs. Rather than looking at the moment of onset, it evaluates the final segmented copy number states observed in terminal clones or cells across focal, medium, and broad scales.

- **E1: Are accumulated complex hcCNAs correct?** Compares the final predicted copy number vectors against the ground truth for terminal cells or clones. It uses RMSE and SCC to evaluate the similarity of the profiles and ACC to measure the exact match rate of the phased CN strings.
- **E2: Are accumulated Mirrored-Subclonal-CNAs correct?** Specifically evaluates Mirrored-Subclonal events—patterns where different subclones have reciprocal alterations—at the clone level to determine if these complex evolutionary outcomes are accurately reconstructed.



```python
# E1
bench.segmentation(
    tool_cna_files=[
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    gt_cna_file="/path/to/gt/cell/haplotype_combined.csv",
    tool_names=[
        "signals",
        "seacon",
    ],
    threshold=0.8,            # Similarity threshold used by the segmentation evaluation
    profile_bin_size=100000,  # Bin size used in the GT profile (e.g., 100kb)
    outprefix="segmentation_thre0.8",
)
# E2
bench..mirrorsubclone(
    tool_cna_files=[
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    tool_names=[
        "signals",
        "seacon",
    ],
    changes_file= "/path/to/gt/mirrored_subclonal_cnas.csv",
    profile_bin_size=100000,
    cal_event_level = True
)

```

------

## Perspective F: Global CN Phasing

This perspective evaluates the accuracy of haplotype-specific phasing, determining how effectively a model assigns copy number states to the correct parental haplotype (e.g., distinguishing a 2|1 state from 1|2).

- **F1: Are CN values correctly phased?** Measures the mismatch error between predicted and ground-truth phase strings. Evaluation is performed in two modes:

  - **GT-heterozygous mode:** Reflects combined errors in both phasing and the ability to retain heterozygous loci.

  - **Shared-heterozygous mode:** Isolates the phasing accuracy specifically at loci where both the prediction and ground truth agree the site is informative (heterozygous).

    The analysis uses normalized Hamming distance to account for potential global phase flips.

```python
# Shared-heterozygous mode
bench.hcPhasing(
    tool_hap1_cna_files=[
        "/path/to/signals/haplotype_1.csv",
        "/path/to/seacon/haplotype_1.csv",
    ],
    tool_hap2_cna_files=[
        "/path/to/signals/haplotype_2.csv",
        "/path/to/seacon/haplotype_2.csv",
    ],
    tool_names=[
        "signals",
        "seacon",
    ],
    ground_truth_hap1_file="/path/to/gt/cell/haplotype_1.csv",
    ground_truth_hap2_file="/path/to/gt/cell/haplotype_2.csv",
    mode="heterozygous-only",
    is_clone=False,
    # outprefix="hcPhasing_both",        # Uncomment if your implementation supports outprefix
)
# GT-heterozygous mode
bench.hcPhasing(
    tool_hap1_cna_files=[
        "/path/to/signals/haplotype_1.csv",
        "/path/to/seacon/haplotype_1.csv",
    ],
    tool_hap2_cna_files=[
        "/path/to/signals/haplotype_2.csv",
        "/path/to/seacon/haplotype_2.csv",
    ],
    tool_names=[
        "signals",
        "seacon",
    ],
    ground_truth_hap1_file="/path/to/gt/cell/haplotype_1.csv",
    ground_truth_hap2_file="/path/to/gt/cell/haplotype_2.csv",
    mode="homozygous-inclusive-pred",
    is_clone=False,
    # outprefix="hcPhasing_gt",
)
```





## Overview of Submodules

| Submodule                                           | Description                                                  | Core Metrics                             |
| --------------------------------------------------- | ------------------------------------------------------------ | ---------------------------------------- |
| [clusterConsistency](clusterConsistency.md)         | Provide a global assessment of clone detection               | AMI, ARI                                 |
| [cloneSizebycluster](cloneSizebycluster.md)         | Evaluate the precision of clone size estimation by comparing predicted cell counts against Ground Truth across different clone scales. | Predicted Size                           |
| [cellprofile](cellprofile.md)                       | Assess the heterogeneity and resolution of CNA calling by counting the number of unique genomic profiles identified. | Unique Group Profile Count               |
| [cloneSizebycellprofile](cloneSizebycellprofile.md) |                                                              |                                          |
| [dolazactree](dolazactree.md)                                         |                                                              |                                          |
| [cndetect](cndetect.md)                                         | Evaluates copy number (CN) detection accuracy                | RMSE, ACC, SCC                           |
| [cnclass](cnclass.md)                                           | Calculates CN state classification metrics                   | AUROC, AUPRC, ACC, Precision, Recall, F1 |
| [hconsetcn](hconsetcn.md)                                         | Evaluates accuracy of inferred parental CN                   | RMSE, ACC                                |
| [hconsetacc](hconsetacc.md)                                     | Checks if CN events are detected at the correct branch of evolution | ACC                                      |
| [hccnchange](hccnchange.md)                                        | Evaluates correctness of CN changes between parent and child clones | RMSE, ACC                                |
| [hccnstable](hccnstable.md)                                        | Calculates evolutionary CN stability along the phylogenetic tree | ACC                                      |
| [segmentation](segmentation.md)                                        |                                                              |                                          |
| [mirrorsubclone](mirrorsubclone.md)                                    | Evaluates mirror-subclone CNA detection accuracy             | RMSE, ACC                                |
| [hcPhasing](hcPhasing.md)                                         | Evaluates haplotype phasing accuracy                         | Mismatch Error, Switch Error             |

------









