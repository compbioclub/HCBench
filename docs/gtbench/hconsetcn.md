## Function

```python
hconsetcn(
    self,
    tool_hap1_cna_files: List[str],
    tool_hap2_cna_files: List[str],
    tool_names: List[str],
    changes_file: str,
    outfile: str = "evolution_onset_parent_CN.csv",
    profile_bin_size=100000
) -> pd.DataFrame
```

This function evaluates **evolution onset (parent → child) copy-number change accuracy** for multiple tools at a specified bin resolution.

Given a ground-truth “change list” (segments with known CNA changes and labels) and each tool’s haplotype-specific CNA profiles, it:

1. Re-bins the ground-truth segments and tool CNA profiles to `profile_bin_size`.
2. Joins the change segments with predicted haplotype CN values (hap1 and hap2).
3. For each change **Type** (e.g., DEL / LOH / gain classes, depending on your `changes_file`), computes:
   - **ACC**: accuracy of parent CN prediction under onset constraints
   - **RMSE**: parent–child CN error (from your internal RMSE routine)
4. Outputs a per-tool, per-type summary table.

------

## Parameters

| Name                  | Type        | Description                                                  |
| --------------------- | ----------- | ------------------------------------------------------------ |
| `tool_hap1_cna_files` | `List[str]` | List of tool CNA profile CSVs for haplotype 1. Must align with `tool_names`. |
| `tool_hap2_cna_files` | `List[str]` | List of tool CNA profile CSVs for haplotype 2. Must align with `tool_names`. |
| `tool_names`          | `List[str]` | Tool names used in result rows.                              |
| `changes_file`        | `str`       | Path to the ground-truth change table (segments + change labels). |
| `outfile`             | `str`       | Output filename for the summary table. Default: `"evolution_onset_parent_CN.csv"`. |
| `profile_bin_size`    | `int`       | Bin size used to split/standardize both GT segments and predictions. Default: `100000` (100kb). |

------

## Input File Format

### `changes_file` (ground-truth change table)

Expected to contain at least:

- `Segment`: genomic interval identifier (used as the region key)
- `Type`: change category label (used for stratified evaluation)
- `Change`: ground-truth change signal used by downstream scoring

Optional:

- `Haplotype`: if present, will be normalized via `self._normalize_hap_label(...)`

Implementation details:

- loaded via `read_and_drop_empty(changes_file)`
- rebinned using:

```python
change_truth_r = split_all_regions(change_truth_r.set_index("Segment"), profile_bin_size)
change_truth_r = change_truth_r.reset_index().rename(columns={"index": "Segment"})
```

So `Segment` must be splittable by your `split_all_regions` logic.

### Tool CNA profiles (`tool_hap1_cna_files`, `tool_hap2_cna_files`)

Expected to be CNA matrices with:

- `region` column
- remaining columns as cell IDs
- values: haplotype-specific copy numbers (format depends on your pipeline)

Each file is:

- loaded via `read_and_drop_empty(...)`
- rebinned and converted back to a `region` column:

```python
p_h1 = split_all_regions(p_h1.set_index("region"), profile_bin_size).reset_index()
p_h2 = split_all_regions(p_h2.set_index("region"), profile_bin_size).reset_index()
```

------

## Evaluation Logic

For each tool and each unique change `Type`:

1. Subset GT changes:

```python
change_truth = change_truth_r[change_truth_r["Type"] == t]
```

1. Join GT change segments with predicted CN for each haplotype:

```python
h1 = self._onset_join(change_truth, p_h1, "hap1")
h2 = self._onset_join(change_truth, p_h2, "hap2")
comb = pd.concat([h1, h2], ignore_index=True)
```

The joined table is expected (by downstream code) to include:

- `Change` (GT)
- `Parent_predict_num` (predicted parent CN)
- `Child_predict_num` (predicted child CN)

1. Compute metrics:

- Accuracy of parent onset prediction:

```python
acc = self._parent_onset_acc(gt, pd_p)
```

- RMSE from parent/child predictions:

```python
pr, _ = self._parent_child_rmse(gt, pd_p, pd_c)
```

1. Append summary row:

- `Tool`
- `Type`
- `RMSE` (parent-related RMSE returned as `pr`)
- `ACC`

------

## Output

A single CSV is written to `self.output_dir`:

```
os.path.join(self.output_dir, outfile)
```

Default:

```
{self.output_dir}/evolution_onset_parent_CN.csv
```

The output table has one row per **(Tool × Type)**.

### Output Columns

| Column | Meaning                                                      |
| ------ | ------------------------------------------------------------ |
| `Tool` | Tool name                                                    |
| `Type` | Change type/category from `changes_file`                     |
| `RMSE` | Parent/parent-child RMSE returned by `self._parent_child_rmse` |
| `ACC`  | Parent onset accuracy returned by `self._parent_onset_acc`   |

------

## Return Value

Returns a `pd.DataFrame` containing the same summary table saved to CSV.

------

## Example

```python
from hcbench.gtbench.gtbench import GTBench

bench = GTBench(output_dir="out/gt_output")

df = bench.hconsetcn(
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
    outfile="evolution_onset_parent_CN.csv",
)

print(df)
```

