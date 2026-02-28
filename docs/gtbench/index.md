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

bench = GTBench(output_dir="/path/to/results")

# 1. CN Detection Benchmark
bench.cndetect(
    tool_cna_files=[
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    cna_profile_file="/path/to/gt/ground_truth_combined_cnv.csv",
    tool_names=["signals", "seacon"]
)

# 2. CN Classification Benchmark
bench.cnclass(
    tool_hap1_cna_files=[
        "/path/to/signals/haplotype_1.csv",
        "/path/to/seacon/haplotype_1.csv",
    ],
    tool_hap2_cna_files=[
        "/path/to/signals/haplotype_2.csv",
        "/path/to/seacon/haplotype_2.csv",
    ],
    tool_names=["signals", "seacon"],
    profile_hap1_cna_file="/path/to/gt/haplotype_1.csv",
    profile_hap2_cna_file="/path/to/gt/haplotype_2.csv",
    type="hcCNA"
)

# 3. CN Change Benchmark
bench.hccnchange(
    tool_hap1_cna_files=[
        "/path/to/signals/haplotype_1.csv",
        "/path/to/seacon/haplotype_1.csv",
    ],
    tool_hap2_cna_files=[
        "/path/to/signals/haplotype_2.csv",
        "/path/to/seacon/haplotype_2.csv",
    ],
    tool_names=["signals", "seacon"],
    changes_file="/path/to/gt/changes.csv",
)

# 4. CN Stability Benchmark
bench.hccnstable(
    tool_hap1_cna_files=[
        "/path/to/signals/haplotype_1.csv",
        "/path/to/seacon/haplotype_1.csv",
    ],
    tool_hap2_cna_files=[
        "/path/to/signals/haplotype_2.csv",
        "/path/to/seacon/haplotype_2.csv",
    ],
    tool_names=["signals", "seacon"],
    changes_file="/path/to/gt/changes.csv",
    tree_file="/path/to/gt/tree.newick",
)

# 5. CN Onset Accuracy Benchmark
bench.hconsetacc(
    tool_hap1_cna_files=[
        "/path/to/signals/haplotype_1.csv",
        "/path/to/seacon/haplotype_1.csv",
    ],
    tool_hap2_cna_files=[
        "/path/to/signals/haplotype_2.csv",
        "/path/to/seacon/haplotype_2.csv",
    ],
    tool_names=["signals", "seacon"],
    changes_file="/path/to/gt/changes.csv",
)

# 6. Parental CN Accuracy Benchmark
bench.hconsetcn(
    tool_hap1_cna_files=[
        "/path/to/signals/haplotype_1.csv",
        "/path/to/seacon/haplotype_1.csv",
    ],
    tool_hap2_cna_files=[
        "/path/to/signals/haplotype_2.csv",
        "/path/to/seacon/haplotype_2.csv",
    ],
    tool_names=["signals", "seacon"],
    changes_file="/path/to/gt/changes.csv",
)

# 7. Haplotype Phasing Accuracy
bench.hcPhasing(
    tool_hap1_cna_files=[
        "/path/to/signals/haplotype_1.csv",
        "/path/to/seacon/haplotype_1.csv",
    ],
    tool_hap2_cna_files=[
        "/path/to/signals/haplotype_2.csv",
        "/path/to/seacon/haplotype_2.csv",
    ],
    tool_names=["signals", "seacon"],
    ground_truth_hap1_file="/path/to/gt/haplotype_1.csv",
    ground_truth_hap2_file="/path/to/gt/haplotype_2.csv",
)

# 8. Mirror-Subclone Benchmark
bench.mirrorsubclone(
    tool_cna_files=[
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    tool_names=["signals", "seacon"],
    changes_file="/path/to/gt/mirrored_clones.csv",
)
```

------

























