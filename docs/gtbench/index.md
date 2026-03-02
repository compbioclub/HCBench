# GTBench

`GTBench` is an integrated benchmarking module for evaluating copy number alteration (CNA) results at both single-cell and clone levels.
 It provides multiple submodules that assess copy number detection, classification, evolutionary stability, haplotype phasing, and clustering performance.

------

## Overview of Submodules

| Submodule                                           | Description                                                  | Core Metrics                             |
| --------------------------------------------------- | ------------------------------------------------------------ | ---------------------------------------- |
| [clusterConsistency](clusterConsistency.md)         | Provide a global assessment of clone detection               | AMI, ARI                                 |
| [cloneSizebycluster](cloneSizebycluster.md)         | Evaluate the precision of clone size estimation by comparing predicted cell counts against Ground Truth across different clone scales. | Predicted Size                           |
| [cellprofile](cellprofile.md)                       | Assess the heterogeneity and resolution of CNA calling by counting the number of unique genomic profiles identified. | Unique Group Profile Count               |
| [cloneSizebycellprofile](cloneSizebycellprofile.md) |                                                              |                                          |
|                                                     |                                                              |                                          |
|                                                     |                                                              |                                          |
|                                                     |                                                              |                                          |
|                                                     |                                                              |                                          |
|                                                     |                                                              |                                          |
| `cndetect`                                          | Evaluates copy number (CN) detection accuracy                | RMSE, ACC, SCC                           |
| `cnclass`                                           | Calculates CN state classification metrics                   | AUROC, AUPRC, ACC, Precision, Recall, F1 |
| `hccnchange`                                        | Evaluates correctness of CN changes between parent and child clones | RMSE, ACC                                |
| `hccnstable`                                        | Calculates evolutionary CN stability along the phylogenetic tree | ACC                                      |
| `hconsetacc`                                        | Checks if CN events are detected at the correct branch of evolution | ACC                                      |
| `hconsetcn`                                         | Evaluates accuracy of inferred parental CN                   | RMSE, ACC                                |
| `hcPhasing`                                         | Evaluates haplotype phasing accuracy                         | Mismatch Error, Switch Error             |
| `mirrorsubclone`                                    | Evaluates mirror-subclone CNA detection accuracy             | RMSE, ACC                                |

------

## Example Usage

Below is a complete example showing how to initialize and run all major GTBench modules.

```
from hcbench.gtbench.gtbench import GTBench

# All benchmark outputs will be written under this directory
bench = GTBench(output_dir="/path/to/results/gt_output/cell_level")

# ============================================================
# 1) cellprofile: Count unique cell-level CNA profiles per tool
# ============================================================
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

# (Optional) Use minor/major (accNA-style) profiles instead of haplotype-combined
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


# ============================================================
# 2) dolazactree: Compute parsimony scores (include GT as baseline)
# ============================================================
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

# (Optional) Parsimony on minor/major (accNA-style) profiles
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


# ============================================================
# 3) clusterConsistency
# ============================================================
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


# ============================================================
# 4) cloneSizebycellprofile: Evaluate clone-size accuracy based on cell profiles
# ============================================================
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

# (Optional) Clone-size accuracy using minor/major (accNA-style) profiles
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


# ============================================================
# 5) cloneSizebycluster: Evaluate clone-size accuracy based on cluster assignments
# ============================================================
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


# ============================================================
# 6) cndetect: Bin-level CN detection benchmark against a CNA profile (GT)
# ============================================================
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

# (Optional) CN detection using minor/major (accNA-style) profiles
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


# ============================================================
# 7) cnclass: CN event classification benchmark (requires hap1 + hap2 inputs)
# ============================================================
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

# (Optional) Treat minor as "hap1" and major as "hap2" for accNA-style evaluation
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


# ============================================================
# 8) hcPhasing: Evaluate haplotype phasing accuracy (two common modes)
# ============================================================
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


# ============================================================
# 9) Evolution-onset benchmarks: run for focal / medium / broad changes
#    - hccnchange: CN change detection
#    - hccnstable: CN stability accuracy (requires a phylogenetic tree)
#    - hconsetacc: onset accuracy
#    - hconsetcn: parent CN correctness at onset
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


# ============================================================
# 10) NA_ratio: Compute the fraction of NA entries in each tool output
# ============================================================
bench.NA_ratio(
    tool_cna_files=[
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    tool_names=[
        "signals",
        "seacon",
    ],
    outfile="NA_ratio.csv",
)


# ============================================================
# 11) segmentation: Evaluate segmentation consistency against GT
# ============================================================
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


# ============================================================
# 12) detect_size: Detection accuracy stratified by event size (cell-level / clone-level)
# ============================================================
bench.detect_size(
    tool_cna_files=[
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    tool_names=[
        "signals",
        "seacon",
    ],
    level="clone",
    size_file="/path/to/gt/clone/ground_truth_classified.csv",
    outfile="detect_size_acc_clone.csv",
)

bench.detect_size(
    tool_cna_files=[
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    tool_names=[
        "signals",
        "seacon",
    ],
    level="cell",
    size_file="/path/to/gt/cell/ground_truth_classified.csv",
    outfile="detect_size_acc_cell.csv",
)


```

------

























