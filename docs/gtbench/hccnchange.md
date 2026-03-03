## Function

```python
hccnchange(
    self,
    tool_hap1_cna_files: List[str],
    tool_hap2_cna_files: List[str],
    tool_names: List[str],
    changes_file: str,
    profile_bin_size=100000,
    outfile: str = "evolution_onset_CN_Change.csv",
) -> pd.DataFrame
```

This function evaluates **haplotype-specific copy-number change (CN change) prediction** for multiple tools during evolutionary onset.

Given a ground-truth change table (`changes_file`) and each tool’s haplotype CNA profiles (hap1 and hap2), it:

1. Standardizes GT segments and tool CNA profiles to the same bin size (`profile_bin_size`).
2. Joins GT change segments with predicted haplotype copy numbers using `_onset_join`.
3. Produces predicted CN-change labels/values (`Change_predict`) and compares them to GT `Change`.
4. Computes **RMSE** and **ACC** for CN-change prediction, stratified by change `Type`.
5. Saves a per-tool, per-type summary table.

------

## Parameters

| Name                  | Type        | Description                                                  |
| --------------------- | ----------- | ------------------------------------------------------------ |
| `tool_hap1_cna_files` | `List[str]` | List of tool CNA profile CSVs for haplotype 1. Must align with `tool_names`. |
| `tool_hap2_cna_files` | `List[str]` | List of tool CNA profile CSVs for haplotype 2. Must align with `tool_names`. |
| `tool_names`          | `List[str]` | Tool names used in output rows.                              |
| `changes_file`        | `str`       | Path to the ground-truth CN-change table (segments + change labels). |
| `profile_bin_size`    | `int`       | Bin size used to split/standardize both GT segments and predictions. Default: `100000` (100kb). |
| `outfile`             | `str`       | Output filename for the summary table. Default: `"evolution_onset_CN_Change.csv"`. |

------

## Input File Format

### `changes_file` (ground-truth CN-change table)

Expected to contain at least:

- `Segment`: genomic interval identifier used as the segment key
- `Type`: change category label (used for stratified evaluation)
- `Change`: ground-truth CN-change signal (numeric or categorical, depending on your pipeline)

Optional:

- `Haplotype`: if present, haplotype labels are normalized via `self._normalize_hap_label(...)`

Implementation details:

```python
change_truth = read_and_drop_empty(changes_file)

if "Haplotype" in change_truth.columns:
    change_truth["Haplotype"] = self._normalize_hap_label(change_truth["Haplotype"])

change_truth = split_all_regions(change_truth.set_index("Segment"), profile_bin_size)
change_truth = change_truth.reset_index().rename(columns={"index": "Segment"})
```

### Tool CNA profiles (`tool_hap1_cna_files`, `tool_hap2_cna_files`)

Expected to be CNA matrices with:

- `region` column
- cell columns
- haplotype-specific CN values

Each file is loaded via `read_and_drop_empty` and re-binned:

```python
p_h1 = split_all_regions(p_h1.set_index("region"), profile_bin_size).reset_index()
p_h2 = split_all_regions(p_h2.set_index("region"), profile_bin_size).reset_index()
```

------

## Evaluation Logic

For each tool:

1. Join GT changes with predicted haplotype CN:

```python
h1 = self._onset_join(change_truth, p_h1, "hap1")
h2 = self._onset_join(change_truth, p_h2, "hap2")
comb = pd.concat([h1, h2], ignore_index=True)
```

The joined table is expected to include (at minimum):

- `Type`
- `Change` (GT)
- `Change_predict` (predicted CN-change)

1. For each change `Type`, compute CN-change metrics:

```python
rmse = self._cn_change_rmse(gt, pd_)
acc  = self._cn_change_acc(gt, pd_)
```

1. Append one summary row per `(Tool × Type)`:

- `Tool`
- `Type`
- `RMSE`
- `ACC`

------

## Output

A single CSV summary table is written to `self.output_dir`:

```
os.path.join(self.output_dir, outfile)
```

Default:

```
{self.output_dir}/evolution_onset_CN_Change.csv
```

------

## Output Table Schema

The output table contains one row per **(Tool × Type)**:

| Column | Meaning                                              |
| ------ | ---------------------------------------------------- |
| `Tool` | Tool name                                            |
| `Type` | Change category label from the GT table              |
| `RMSE` | CN-change RMSE computed by `self._cn_change_rmse`    |
| `ACC`  | CN-change accuracy computed by `self._cn_change_acc` |

Notes:

- The exact interpretation of `RMSE` and `ACC` depends on how `Change` / `Change_predict` are encoded and how `_cn_change_rmse` / `_cn_change_acc` are implemented (numeric vs categorical).

------

## Return Value

Returns a `pd.DataFrame` with columns:

- `Tool`
- `Type`
- `RMSE`
- `ACC`

------

## Example

```python
from hcbench.gtbench.gtbench import GTBench

bench = GTBench(output_dir="out/gt_output")

df = bench.hccnchange(
    tool_hap1_cna_files=[
        "/path/to/chisel/hap1.csv",
        "/path/to/signals/hap1.csv",
    ],
    tool_hap2_cna_files=[
        "/path/to/chisel/hap2.csv",
        "/path/to/signals/hap2.csv",
    ],
    tool_names=["CHISEL", "SIGNALS"],
    changes_file="/path/to/gt/changes_cn_change.csv",
    profile_bin_size=100000,
    outfile="evolution_onset_CN_Change.csv",
)

print(df)
```

