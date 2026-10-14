# HCBench-sim

`HCBench-sim` evaluates inferred hcCNA profiles against the hcCNA-bench GT. The evaluation is organized around five evaluation questions: CNA detection, clone detection, evolution tracking, accumulated complex hcCNA detection, and CN phasing.

## Prepare inputs and initialize

Prepare standardized caller CNA matrices and matching HCDSIM ground truth. Combined files use a first `region` column and one column per cell; split files contain one haplotype or allele per matrix. Clone files need `cell_id` and `clone_id`. Evolution-onset methods also need the corresponding GT change tables; stability analysis needs a Newick tree.

Keep caller paths and `tool_names` in the same order. Use cell-level files together or clone-level files together. Replace all `/path/to/` paths below with your data. Run the setup block once before the method examples.

`bin_size` and `genome` are required when creating a CNA parser, `GTBench`, or `RealBench`:

- `bin_size`: a positive integer in base pairs, such as `100000` for 100 kb. Use the resolution of your CNA data.
- `genome`: `"hg19"` or `"hg38"`, matching your input coordinates.

Methods use this configuration. Initialization prints:

```text
Bin size: 100000 bp
Genome: hg38
```

With the default event-size thresholds, analysis is available for:

| CNA `bin_size` | Event categories |
| --- | --- |
| `<= 3,000,000` bp | focal, medium, broad |
| `> 3,000,000` and `<= 10,000,000` bp | medium, broad |
| `> 10,000,000` bp | broad |

Event length is `End - Start + 1`: focal is below 3 Mb, medium is from 3 Mb through 10 Mb, and broad is above 10 Mb. Whole-chromosome and whole-genome events are reported as broad. Categories unavailable at your resolution have no annotation or category-metric rows; a missing comparison is reported as NA rather than evidence of no events.

For evolution-onset analyses, supply change tables containing only the event categories appropriate for your resolution. Those methods read the change tables you provide. Phylogeny examples require `lazac` on `PATH`.

```python
from hcbench.gtbench import GTBench

# All benchmark outputs will be written under this directory
bench = GTBench(output_dir="/path/to/results/gt_output/cell_level", bin_size=100000, genome="hg38")
```

------

All files are saved under `bench.output_dir`. Methods with `outfile` use that filename; phasing uses `{outprefix}_{mode}.csv`. `cna_values()` returns its metrics DataFrame as well as saving the CSV. Methods that only save files, including `complex_hccnas()` and unique-profile counts/sizes, return `None`.

## CNA Detection

This section evaluates the accuracy of allele-specific copy number alterations (acCNAs) and haplotype-specific copy number alterations (hcCNAs) at the value and state levels.

### acCNA values

**Are allele-, cell/clone-, and locus-specific CN values correct?**

Assesses quantitative accuracy using RMSE and SCC for major and minor alleles. It also uses ACC to measure the percentage of loci where the full allele-specific string is exactly correct.

```python
# acCNA values
bench.cna_values(
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
    outfile="cna_values_accna.csv",
)
```

### acCNA states

**Are allele-, cell/clone-, and locus-specific CN loss, gain, and neutral states correct?**

Treats copy number detection as a binary classification problem (Gain, Neutral, Loss). Performance is measured via a suite of statistical metrics (AUROC, Sensitivity, F1-Score, etc.) to determine the reliability of calling specific genomic states.

```python
# acCNA states
bench.cna_states(
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
    type="acCNA",
    outfile="cna_states_accna.csv",
)
```

### hcCNA values

**Are haplotype-, cell/clone-, and locus-specific CN values correct?**

Calculates the average RMSE and SCC across both individual haplotypes and evaluates the rate of exact matches for the phased CN strings (ACC).

```python
# hcCNA values
bench.cna_values(
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
    outfile="cna_values_hccna.csv",
)
```

### hcCNA states

**Are haplotype-, cell/clone-, and locus-specific CN loss, gain, and neutral states correct?**

Assesses phased CN calling performance via binary classification across multiple states (Gain/Loss/Neutral), providing a comprehensive statistical breakdown of detection sensitivity and specificity.

```python
# hcCNA states
bench.cna_states(
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
    type="hcCNA",
    outfile="cna_states_hccna.csv",
)
```

------

## Clone Detection

This section evaluates the accuracy of inferring clonal architecture and the diversity of unique CN profiles within a tumor population.

### Global tumor clones

**Are global tumor clones correct?**

Uses global metrics like Adjusted Mutual Information (AMI) and Adjusted Rand Index (ARI) to compare predicted labels against the ground truth. This provides a high-level assessment of how well the overall population structure is captured.

```python
# Global tumor clones
bench.global_tumor_clones(
    tool_clone_files=[
        "/path/to/signals/clusters.csv",
        "/path/to/seacon/clusters.csv",
    ],
    tool_names=[
        "signals",
        "seacon",
    ],
)
```

### Rare and dominant clones

**Are rare and dominant clones correct?**

Evaluates clustering accuracy stratified by clone size using Clone Size Deviation (CSD). It measures whether specific clones (regardless of size) are being merged into larger clusters (positive deviation) or "shattered" into many small groups (negative deviation).

First, generate the GT `clusters.csv` from HCDSIM's `profile/barcodes.txt`. The barcode file has no header and contains one cell ID per line. In HCDSIM cell IDs, the prefix before the first underscore identifies the clone; for example, `clone1_cell1` belongs to `clone1`.

```python
from pathlib import Path
import pandas as pd

df = pd.read_csv(
    "/path/to/HCDSIM/profile/barcodes.txt",
    sep="\t",
    comment=None,
    index_col=False,
    header=None,
    names=["cell_id"],
)
df["clone_id"] = df["cell_id"].str.split("_", n=1).str[0]

gt_cluster_file = Path("results/gt/clusters.csv")
gt_cluster_file.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(gt_cluster_file, index=False)
```

This writes a CSV with `cell_id` and `clone_id` columns. Keep the full cell IDs so they match those in the caller cluster files. Use the generated GT file with the `bench` initialized above:

```python
# Rare and dominant clones
bench.rare_and_dominant_clones(
    tool_cluster_files=[
        "/path/to/signals/clusters.csv",
        "/path/to/seacon/clusters.csv",
    ],
    gt_cluster_file=str(gt_cluster_file),
    tool_names=[
        "signals",
        "seacon",
    ],
    outfile="rare_and_dominant_clones.csv",
)
```

### Unique acCNA profiles

**Are unique allele- and cell-specific CN profiles correct?**

Analyzes "Unique Profile Groups" (UPGs) using Unique Profile Count Error (UPCE) and Unique Profile Size Deviation (UPSD). These metrics detect if distinct biological profiles are being over-smoothed (merged) or over-fragmented due to noise.

```python
# Unique acCNA profiles
bench.unique_cna_profile_counts(
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
    outfile="unique_cna_profile_counts_accna.csv",
)

bench.unique_cna_profile_sizes(
    tool_cna_files=[
        "/path/to/signals/minor_major.csv",
        "/path/to/seacon/minor_major.csv",
    ],
    gt_cna_file="/path/to/gt/cell/minor_major.csv",
    tool_names=[
        "signals",
        "seacon",
    ],
    outfile="unique_cna_profile_sizes_accna.csv",
)
```

### Unique hcCNA profiles

**Are unique haplotype- and cell-specific CN profiles correct?**

Quantifies performance based on cells sharing identical haplotype-specific CN strings across the genome using UPCE and UPSD metrics.

```python
# Unique hcCNA profiles
bench.unique_cna_profile_counts(
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
    outfile="unique_cna_profile_counts_hccna.csv",  # Output CSV name under output_dir
)

bench.unique_cna_profile_sizes(
    tool_cna_files=[
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    gt_cna_file="/path/to/gt/cell/haplotype_combined.csv",
    tool_names=[
        "signals",
        "seacon",
    ],
    outfile="unique_cna_profile_sizes_hccna.csv",
)
```

------

## Evolution Tracking

This section evaluates whether inferred CN phylogenies and complex hcCNA evolutionary events agree with the ground truth.

### acCNA phylogeny

**Are allele- and cell-specific CN phylogenies correct?**

Calculates a parsimony score for each supplied CN profile. Include GT in the input list to compare its score with caller scores.

```python
# acCNA phylogeny
bench.cna_phylogeny(
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
    outfile="cna_phylogeny_accna.csv",
)
```

### hcCNA phylogeny

**Are haplotype- and cell-specific CN phylogenies correct?**

Calculates parsimony scores from phased CN profiles. The output lists each tool and its `Parsimony` score.

```python
# hcCNA phylogeny
bench.cna_phylogeny(
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
    outfile="cna_phylogeny_hccna.csv",
)
```

### Prepare size-specific change tables

Before running the four onset and stability analyses below, split the HCDSIM `cell_changes.csv` and `clone_changes.csv` tables by event size. Use `Length` in base pairs: focal is below 3,000,000 bp, medium is from 3,000,000 through 10,000,000 bp, and broad is above 10,000,000 bp. Keep the other change-table columns for the benchmark methods.

The following example reads both files from your HCDSIM profile directory and saves one CSV per size present in each table. Replace the input directory with yours. Missing or nonnumeric lengths must be corrected before classification.

```python
from pathlib import Path
import numpy as np
import pandas as pd

profile_dir = Path("/path/to/HCDSIM/profile")
changes_by_level = {}

for level in ("cell", "clone"):
    df = pd.read_csv(profile_dir / f"{level}_changes.csv")
    df["Length"] = pd.to_numeric(df["Length"], errors="coerce")
    if df["Length"].isna().any():
        raise ValueError(f"Correct missing or nonnumeric Length values in {level}_changes.csv")

    df["size"] = np.select(
        [
            df["Length"] < 3_000_000,
            (df["Length"] >= 3_000_000) & (df["Length"] <= 10_000_000),
        ],
        ["focal", "medium"],
        default="broad",
    )

    output_dir = Path("results/gt") / f"{level}-level"
    output_dir.mkdir(parents=True, exist_ok=True)
    changes_by_level[level] = {}
    if level == "cell":
        df.rename(columns = {
            "Clone": "Parent",
            "Cell": "Child"
        }, inplace=True)

    for size in df["size"].unique():
        out_file = output_dir / f"{level}_changes_split_annotated_{size}.csv"
        df.loc[df["size"] == size].to_csv(out_file, index=False)
        changes_by_level[level][size] = out_file
```

Cell-level files are saved under `results/gt/cell-level/` with names such as `cell_changes_split_annotated_focal.csv`; clone-level files are saved under `results/gt/clone-level/` with names such as `clone_changes_split_annotated_focal.csv`. Categories absent from an input table produce no file.

The examples below use the cell-level files. For clone-level analysis, select `changes_by_level["clone"]` and supply matching clone-level haplotype matrices. For stability analysis, the tree node names must match the parent and child identifiers in the chosen change table. These methods evaluate the supplied change file; use one size-specific file per call.

### Parental CN at onset

**Are parent CNs of complex hcCNAs correct at onset?** Assesses whether the model correctly identifies the pre-event (parental) copy number state at the exact moment an alteration arises during evolution, using metrics like ACC and RMSE.

### CN change at onset

**Are CN changes of complex hcCNAs correct at onset?** Evaluates the model's accuracy in capturing the magnitude and direction of the copy number transition (e.g., a gain from 1 to 3) from parent to child.

### Complex hcCNA at onset

**Are complex hcCNAs correct at onset?** A holistic measure requiring the parent state, child state, and the resulting change to all match the ground truth. This reflects the end-to-end reliability in pinpointing true evolutionary transitions.

### Complex hcCNA stability

**Are complex hcCNAs stable throughout evolution?** Evaluates whether a detected alteration is correctly tracked as a stable, heritable state across all descendant clones following its initial onset.

Run the preparation step above first. This loop evaluates only sizes that have a change file and are supported by `bench.bin_size`: all three sizes up to 3 Mb, medium/broad above 3 Mb through 10 Mb, and broad above 10 Mb.

```python
# Run after preparing the size-specific change tables above
# Parental CN at onset: parental_cn_at_onset
# CN change at onset: cn_change_at_onset
# Complex hcCNA at onset: complex_hccna_at_onset
# Complex hcCNA stability: complex_hccna_stability
change_files = changes_by_level["cell"]
event_sizes = ["broad"]
if bench.bin_size <= 10_000_000:
    event_sizes.insert(0, "medium")
if bench.bin_size <= 3_000_000:
    event_sizes.insert(0, "focal")
event_sizes = [size for size in event_sizes if size in change_files]
for size in event_sizes:

    bench.parental_cn_at_onset(
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
        changes_file=str(change_files[size]),
        outfile=f"parental_cn_at_onset_{size}.csv",
    )

    bench.cn_change_at_onset(
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
        changes_file=str(change_files[size]),
        outfile=f"cn_change_at_onset_{size}.csv",
    )

    bench.complex_hccna_at_onset(
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
        changes_file=str(change_files[size]),
        outfile=f"complex_hccna_at_onset_{size}.csv",
    )

    bench.complex_hccna_stability(
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
        changes_file=str(change_files[size]),
        tree_file="/path/to/tree/tree_new.newick",
        outfile=f"complex_hccna_stability_{size}.csv",
    )
```

------

## Accumulated Complex hcCNA Detection

This section assesses the final "burden" of complex hcCNAs. Rather than looking at the moment of onset, it evaluates the final segmented copy number states observed in terminal clones or cells across focal, medium, and broad scales.

### Complex hcCNAs

**Are accumulated complex hcCNAs correct?**

Compares predicted complex hcCNA profiles against the ground truth (GT) using RMSE, SCC, and ACC to assess numerical agreement and exact matches of haplotype-specific CN values.

Use the `bench` initialized above. Supply either a combined GT CNA matrix through `gt_cna_file`, or an existing GT annotation DataFrame through `gt_annotated_df`. Keep `tool_cna_files` and `tool_names` in the same order.

#### Evaluate from a GT CNA matrix

When you have a combined GT matrix, `complex_hccnas()` annotates it before comparing the caller profiles:

```python
bench.complex_hccnas(
    tool_cna_files=[
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    gt_cna_file="/path/to/gt/cell/haplotype_combined.csv",
    tool_names=["signals", "seacon"],
    threshold=0.8,
    outprefix="complex_hccnas_",
)
```

The default `threshold=0.8` controls whole-chromosome/genome event annotation. Use `chromosomes=["chr1", "chr2"]` to restrict GT annotation and caller profiles to selected chromosomes.

#### Evaluate with existing GT annotations

If you already have a GT annotation CSV, load it with pandas and pass the resulting DataFrame as `gt_annotated_df`. Set `gt_cna_file=None`; this argument must still be supplied. The method uses the supplied event labels directly, so `threshold` does not change them.

The table needs `cell`, `value`, `size`, and `type`, plus either `region` or all three coordinate columns `Chrom`, `Start`, and `End`. Values use haplotype pairs such as `2|1`; size labels are `focal`, `medium`, or `broad`. Each row describes a cell and an input CNA region, with cell IDs and region coordinates matching the caller matrices at the configured bin resolution. A `Segment` column alone cannot replace these region coordinates.

To create a suitable CSV from a combined GT matrix, save annotations with `detect_cnv_type=True`:

```python
from pathlib import Path
from hcbench.parsers.gt import GTParser

Path("results/gt").mkdir(parents=True, exist_ok=True)
gt_parser = GTParser(
    input_path="/path/to/gt/cell/haplotype_combined.csv",
    output_path="results/gt",
    bin_size=bench.bin_size,
    genome=bench.genome,
)
gt_parser.annotate_segments(
    "/path/to/gt/cell/haplotype_combined.csv",
    threshold=0.8,
    detect_cnv_type=True,
    out_file="ground_truth_classified.csv",
)
```

This writes `results/gt/ground_truth_classified.csv`. `GTParser.annotate_segments()` returns `None`, so load the saved CSV for evaluation:

```python
import pandas as pd

gt_annotated_df = pd.read_csv(
    "results/gt/ground_truth_classified.csv",
    dtype={"cell": str, "Chrom": str, "value": str},
)
bench.complex_hccnas(
    gt_cna_file=None,
    gt_annotated_df=gt_annotated_df,
    tool_cna_files=[
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    tool_names=["signals", "seacon"],
    outprefix="complex_hccnas_",
)
```

When supplying `gt_annotated_df` together with `chromosomes`, prepare the DataFrame with only those chromosomes; the method applies `chromosomes` to the caller profiles but does not filter the supplied GT table by chromosome. 

Both input options write `{outprefix}metrics.csv` under `bench.output_dir` and return `None`. With the examples above, the filename is `complex_hccnas_metrics.csv`. Rows report `Size`, `Type`, `ACC`, `SCC`, `RMSE`, and `Tool`. Event categories excluded by `bench.bin_size` have no metric rows; comparisons without matching cell/region entries have missing metrics.


### Mirrored-subclonal-CNAs

**Are accumulated Mirrored-Subclonal-CNAs correct?** Specifically evaluates Mirrored-Subclonal events—patterns where different subclones have reciprocal alterations—at the clone level to determine if these complex evolutionary outcomes are accurately reconstructed.

```python
# Mirrored-subclonal-CNAs
bench.mirrored_subclonal_cnas(
    tool_cna_files=[
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    tool_names=[
        "signals",
        "seacon",
    ],
    changes_file= "/path/to/gt/mirrored_subclonal_cnas_cnas.csv"
)
```

------

## CN Phasing

This section evaluates the accuracy of haplotype-specific phasing, determining how effectively a model assigns copy number states to the correct parental haplotype (e.g., distinguishing a 2|1 state from 1|2).

### Global CN phasing

**Are CN values correctly phased?** Measures the mismatch error between predicted and ground-truth phase strings. Evaluation is performed in two modes:

- **GT-heterozygous mode:** Reflects combined errors in both phasing and the ability to retain heterozygous loci.
- **Shared-heterozygous mode:** Isolates the phasing accuracy specifically at loci where both the prediction and ground truth agree the site is informative (heterozygous).

The analysis uses normalized Hamming distance to account for potential global phase flips.

```python
# Shared-heterozygous mode
bench.global_cn_phasing(
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
    mode="Shared-heterozygous",
    is_clone=False,
    outprefix="global_cn_phasing_both",
)

# GT-heterozygous mode
bench.global_cn_phasing(
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
    mode="GT-heterozygous",
    is_clone=False,
    outprefix="global_cn_phasing_gt",
)
```
