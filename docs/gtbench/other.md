



## Submodule Descriptions

### 1. cna_values

Calculates **RMSE**, **ACC**, and **SCC** for evaluating CN detection at the bin level.
 It compares predicted CN profiles from different tools to the ground truth.

**Output:** `cna_values.csv`


------

### 2. cna_states

Evaluates the accuracy of CN state classification across bins.
 Supports both `acCNA` and `hcCNA` modes and computes AUROC, AUPRC, accuracy, precision, recall, and F1 score for each tool.

**Output:** Categorized CNV CSV files and evaluation summary tables.

------

### 3. cn_change_at_onset

Determines whether CN changes between parent and child clones are correctly identified.
 Evaluates the **Root Mean Squared Error (RMSE)** and **Accuracy (ACC)** for each tool.

**Output:** `cn_change_at_onset.csv`

------

### 4. complex_hccna_stability

Computes **Evolutionary CN Stability** based on phylogenetic tree relationships between clones.
 Requires both maternal/paternal CNA files and a `tree.newick` file describing clone hierarchy.

**Output:** `complex_hccna_stability.csv`


------

### 5. hconsetccc

Determines whether CN events are correctly detected at their **first evolutionary branch**.
 This submodule evaluates **accuracy (ACC)** across DEL/DUP events.

**Output:** `complex_hccna_at_onset.csv`


------

### 6. parental_cn_at_onset

Assesses whether **parental CN states** are correctly inferred during evolution.
 Evaluates both **RMSE** and **ACC**.

**Output:** `parental_cn_at_onset.csv`


------

### 7. hcphasing

Quantifies **haplotype phasing accuracy** using two metrics:

- **Mismatch Error**: Fraction of individual loci where predicted haplotypes differ from the ground truth.
- **Switch Error**: Frequency of phase switches between consecutive loci.

**Output:** `global_cn_phasing.csv`


------

### 8. mirrored_subclonal_cnas

Evaluates accuracy of detecting **mirror-subclone CNAs**, using RMSE and ACC.
 Mirror-subclones are pairs of subclones that exhibit complementary CNA patterns.

**Output:** `mirrored_subclonal_cnas_cnas.csv`

------

### 9. subdetect

Calculates **Adjusted Mutual Information (AMI)** and **Adjusted Rand Index (ARI)** to measure clustering accuracy of predicted clone structures.
 Each classification file should contain `cell_id` and `clone_id` columns.

**Output:** `AMI_ARI_results.csv`

------

