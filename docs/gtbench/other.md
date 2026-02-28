



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

