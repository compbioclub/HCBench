## Function

```python
hconsetacc(
    self,
    tool_hap1_cna_files: List[str],
    tool_hap2_cna_files: List[str],
    tool_names: List[str],
    changes_file: str,
    profile_bin_size=100000,
    outfile: str = "evolution_onset_acc.csv",
) -> pd.DataFrame
```

This function evaluates **evolution onset change classification accuracy** for multiple tools.

Given a ground-truth change table (segments annotated with change `Type` and GT `Change`) and each tool’s haplotype-specific CNA profiles, it:

1. Re-bins GT change segments and tool CNA profiles to a uniform resolution (`profile_bin_size`).
2. Joins GT segments with tool predictions for **hap1** and **hap2**.
3. Derives a predicted change label (`Change_predict`) via your internal onset-join logic.
4. Computes **accuracy (ACC)** of `Change_predict` vs GT `Change`, stratified by change `Type`.
5. Saves both intermediate combined tables (per tool) and the final summary table.

------

## Parameters

| Name                  | Type        | Description                                                  |
| --------------------- | ----------- | ------------------------------------------------------------ |
| `tool_hap1_cna_files` | `List[str]` | List of tool CNA profile CSVs for haplotype 1. Must align with `tool_names`. |
| `tool_hap2_cna_files` | `List[str]` | List of tool CNA profile CSVs for haplotype 2. Must align with `tool_names`. |
| `tool_names`          | `List[str]` | Tool names used for output file naming and result rows.      |
| `changes_file`        | `str`       | Path to the ground-truth change table (segments + labels).   |
| `profile_bin_size`    | `int`       | Bin size used to split/standardize both GT segments and predictions. Default: `100000` (100kb). |
| `outfile`             | `str`       | Output filename for the summary table. Default: `"evolution_onset_acc.csv"`. |

------

## Input File Format

### `changes_file` (ground-truth change table)

Expected to contain at least:

- `Segment`: genomic segment identifier (used as join key after re-binning)
- `Type`: change category label (used for stratified reporting)
- `Change`: ground-truth change label used for accuracy evaluation

Optional:

- `Haplotype`: if present, will be normalized via `self._normalize_hap_label(...)`

Implementation details:

```python
change_truth = read_and_drop_empty(changes_file)
change_truth = split_all_regions(change_truth.set_index("Segment"), profile_bin_size)
change_truth = change_truth.reset_index().rename(columns={"index": "Segment"})
```

### Tool CNA profiles

Each file is expected to have:

- `region` column
- remaining columns as cell IDs
- values: haplotype-specific CNA states (format depends on your pipeline)

They are loaded and re-binned as:

```python
p_h1 = split_all_regions(p_h1.set_index("region"), profile_bin_size).reset_index()
p_h2 = split_all_regions(p_h2.set_index("region"), profile_bin_size).reset_index()
```

------

## Evaluation Logic

For each tool:

1. Join GT changes with haplotype predictions:

```python
h1 = self._onset_join(change_truth, p_h1, "hap1")
h2 = self._onset_join(change_truth, p_h2, "hap2")
comb = pd.concat([h1, h2], ignore_index=True)
```

The joined table is expected to include (at minimum):

- `Type`
- `Change` (GT)
- `Change_predict` (predicted change label)

1. Save the per-tool combined table for debugging/inspection:

```
{self.output_dir}/{ToolName}_comb_combined.csv
```

1. Compute accuracy per `Type`:

- Compare `Change_predict` vs `Change` only where both are non-missing:

```python
mask = gt.notna() & pd_.notna()
acc = (pd_[mask] == gt[mask]).mean()
```

1. Append one result row per `(Tool × Type)`:

- `Tool`
- `Type`
- `ACC`

------

## Output

### Intermediate per-tool combined tables

For each tool, the function writes:

```
os.path.join(self.output_dir, f"{name}_comb_combined.csv")
```

These files contain the merged GT + predictions across both haplotypes.

### Final summary table (`outfile`)

Saved to:

```
os.path.join(self.output_dir, outfile)
```

Default:

```
{self.output_dir}/evolution_onset_acc.csv
```

------

## Output Table Schema

The final CSV / dataframe contains one row per **(Tool × Type)**:

| Column | Meaning                                                      |
| ------ | ------------------------------------------------------------ |
| `Tool` | Tool name                                                    |
| `Type` | Change category label from the GT table                      |
| `ACC`  | Accuracy of predicted change labels (`Change_predict`) vs GT (`Change`) |

Notes:

- If a given `Type` has no valid comparable entries (all missing), `ACC` is set to `NaN`.

------

## Return Value

Returns a `pd.DataFrame` with columns:

- `Tool`
- `Type`
- `ACC`

------

## Example

```python
from hcbench.gtbench.gtbench import GTBench

bench = GTBench(output_dir="out/gt_output")

df = bench.hconsetacc(
    tool_hap1_cna_files=[
        "/path/to/chisel/hap1.csv",
        "/path/to/signals/hap1.csv",
    ],
    tool_hap2_cna_files=[
        "/path/to/chisel/hap2.csv",
        "/path/to/signals/hap2.csv",
    ],
    tool_names=["CHISEL", "SIGNALS"],
    changes_file="/path/to/gt/changes_onset.csv",
    profile_bin_size=100000,
    outfile="evolution_onset_acc.csv",
)

print(df)
```