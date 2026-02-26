# GTBench

`GTBench` is an integrated benchmarking module for evaluating copy number alteration (CNA) results at both single-cell and clone levels.
 It provides multiple submodules that assess copy number detection, classification, evolutionary stability, haplotype phasing, and clustering performance.

------

## Overview of Submodules

| Submodule        | Description                                                  | Core Metrics                             |
| ---------------- | ------------------------------------------------------------ | ---------------------------------------- |
| [clusterConsistency](clusterConsistency.md) | Provide a global assessment of clone detection | AMI, ARI                    |
| `cndetect`       | Evaluates copy number (CN) detection accuracy                | RMSE, ACC, SCC                           |
| `cnclass`        | Calculates CN state classification metrics                   | AUROC, AUPRC, ACC, Precision, Recall, F1 |
| `hccnchange`     | Evaluates correctness of CN changes between parent and child clones | RMSE, ACC                                |
| `hccnstable`     | Calculates evolutionary CN stability along the phylogenetic tree | ACC                                      |
| `hconsetacc`     | Checks if CN events are detected at the correct branch of evolution | ACC                                      |
| `hconsetcn`      | Evaluates accuracy of inferred parental CN                   | RMSE, ACC                                |
| `hcPhasing`      | Evaluates haplotype phasing accuracy                         | Mismatch Error, Switch Error             |
| `mirrorsubclone` | Evaluates mirror-subclone CNA detection accuracy             | RMSE, ACC                                |

------

## Example Usage

```
from hcbench.gtbench.gtbench import GTBench

bench = GTBench(output_dir="/path/to/output")

bench.cndetect(
    tool_cna_files=["/path/to/signals/haplotype_combined.csv",
                    "/path/to/seacon/haplotype_combined.csv"],
    cna_profile_file="/path/to/ground_truth_combined_cnv.csv",
    tool_names=["signals", "seacon"]
)
```

------

## Submodule Descriptions

### 1. cndetect

Calculates **RMSE**, **ACC**, and **SCC** for evaluating CN detection at the bin level.
 It compares predicted CN profiles from different tools to the ground truth.

**Output:** `bin_level_results.csv`


------

### 2. cnclass

Evaluates the accuracy of CN state classification across bins.
 Supports both `acCNA` and `hcCNA` modes and computes AUROC, AUPRC, accuracy, precision, recall, and F1 score for each tool.

**Output:** Categorized CNV CSV files and evaluation summary tables.

------

### 3. hccnchange

Determines whether CN changes between parent and child clones are correctly identified.
 Evaluates the **Root Mean Squared Error (RMSE)** and **Accuracy (ACC)** for each tool.

**Output:** `evolution_onset_CN_Change.csv`

------

### 4. hccnstable

Computes **Evolutionary CN Stability** based on phylogenetic tree relationships between clones.
 Requires both maternal/paternal CNA files and a `tree.newick` file describing clone hierarchy.

**Output:** `evolution_cn_stability_acc.csv`


------

### 5. hconsetccc

Determines whether CN events are correctly detected at their **first evolutionary branch**.
 This submodule evaluates **accuracy (ACC)** across DEL/DUP events.

**Output:** `evolution_onset_acc.csv`


------

### 6. hconsetcn

Assesses whether **parental CN states** are correctly inferred during evolution.
 Evaluates both **RMSE** and **ACC**.

**Output:** `evolution_onset_parent_CN.csv`


------

### 7. hcphasing

Quantifies **haplotype phasing accuracy** using two metrics:

- **Mismatch Error**: Fraction of individual loci where predicted haplotypes differ from the ground truth.
- **Switch Error**: Frequency of phase switches between consecutive loci.

**Output:** `hcPhasing.csv`


------

### 8. mirrorsubclone

Evaluates accuracy of detecting **mirror-subclone CNAs**, using RMSE and ACC.
 Mirror-subclones are pairs of subclones that exhibit complementary CNA patterns.

**Output:** `mirror_subclone_result.csv`

------

### 9. subdetect

Calculates **Adjusted Mutual Information (AMI)** and **Adjusted Rand Index (ARI)** to measure clustering accuracy of predicted clone structures.
 Each classification file should contain `cell_id` and `clone_id` columns.

**Output:** `AMI_ARI_results.csv`

------



## Running GTBench

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

