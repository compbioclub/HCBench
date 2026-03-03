## Function

```python
mirrorsubclone(
    self,
    tool_cna_files: List[str],
    tool_names: List[str],
    changes_file: str,
    profile_bin_size=100000,
    outfile: str = "mirror_subclone_result.csv",
) -> pd.DataFrame
```

This function evaluates tool performance on **Mirrored-Subclonal CNA (MSCNA)** events by comparing predicted haplotype CNAs for two paired clones against a ground-truth MSCNA table.

Given:

- a ground-truth MSCNA event table (`changes_file`) describing two related subclones (Clone1/Clone2) and their haplotype CN states,
- and each tool’s inferred CNA matrix,

it:

1. Converts MSCNA event intervals into a `region` key and re-bins to `profile_bin_size`.
2. Re-bins tool predictions to the same resolution.
3. Merges GT MSCNA events with tool predictions using `_mirror_merge`.
4. Computes:
   - **ACC**: fraction of events where *all four haplotype CN states* (Clone1 hap1/hap2 + Clone2 hap1/hap2) match exactly.
   - **RMSE**: mean per-event haplotype error aggregated across the four haplotype states.

Results are summarized per tool and written to disk.

------

## Parameters

| Name               | Type        | Description                                                  |
| ------------------ | ----------- | ------------------------------------------------------------ |
| `tool_cna_files`   | `List[str]` | List of tool CNA profile CSV files. Must align with `tool_names`. |
| `tool_names`       | `List[str]` | Tool names used in result rows.                              |
| `changes_file`     | `str`       | Path to the mirrored-subclone ground-truth table (MSCNA events). |
| `profile_bin_size` | `int`       | Bin size used to split/standardize both GT events and predictions. Default: `100000` (100kb). |
| `outfile`          | `str`       | Output filename for the summary table. Default: `"mirror_subclone_result.csv"`. |

------

## Input File Format

### `changes_file` (MSCNA ground-truth table)

Loaded via:

```python
change_df_r = read_and_drop_empty(changes_file)
```

Expected columns include at least:

- `Chromosome`
- `Start`
- `End`

These are combined into a genomic interval string:

```python
change_df_r["region"] = f"{Chromosome}:{Start}-{End}"
```

Additional columns must exist that are required by `_mirror_merge(...)` and subsequent computations.
From the downstream code, the merged table is expected to contain (after `_mirror_merge`) fields like:

- `Clone1_CNA`, `Clone2_CNA` (GT haplotype CNA encoded as `"hap1|hap2"`)
- `Clone1_predict_CNA`, `Clone2_predict_CNA` (predicted haplotype CNA encoded as `"hap1|hap2"`)

GT events are re-binned:

```python
change_truth_r = split_all_regions(change_df_r.set_index("region"), profile_bin_size)
change_truth_r = change_truth_r.reset_index().rename(columns={"index": "region"})
```

### Tool CNA files (`tool_cna_files`)

Expected to be CNA matrices with:

- `region` column
- cell/clone columns depending on your pipeline
- CNA values encoded as haplotype strings like `"a|b"` (required for hap split)

Tool predictions are re-binned similarly:

```python
pred = split_all_regions(pred.set_index("region"), profile_bin_size).reset_index()
```

------

## Processing Steps

For each tool:

1. Merge GT MSCNA table with predictions:

```python
combined = self._mirror_merge(change_df, pred)
```

If `combined` is empty, RMSE and ACC are set to `NaN`.

1. Split haplotype CNA strings for both clones:

For each side in `["Clone1", "Clone2"]`, the function derives:

- `{side}_hap1_CNA`, `{side}_hap2_CNA` (GT)
- `{side}_predict_hap1_CNA`, `{side}_predict_hap2_CNA` (prediction)

by splitting `"hap1|hap2"`.

------

## Metrics

### 1) ACC (Exact match accuracy)

An event is counted as correct only if **all four** haplotype values match exactly:

- Clone1 hap1 and hap2
- Clone2 hap1 and hap2

```python
acc = mean(
    Clone1_pred_h1 == Clone1_gt_h1 AND
    Clone1_pred_h2 == Clone1_gt_h2 AND
    Clone2_pred_h1 == Clone2_gt_h1 AND
    Clone2_pred_h2 == Clone2_gt_h2
)
```

So ACC is the fraction of MSCNA events fully matched across both clones and both haplotypes.

### 2) RMSE (Aggregated per-event haplotype error)

Per haplotype value, squared error is computed after coercing to numeric:

```python
se = (pred - gt)^2
```

Then per-event RMSE-like aggregate is computed as:

```python
RMSE_result =
    sqrt(Clone1_hap1_error) +
    sqrt(Clone1_hap2_error) +
    sqrt(Clone2_hap1_error) +
    sqrt(Clone2_hap2_error)
```

Finally, the reported RMSE is the mean of `RMSE_result` across events.

Note: This is not a conventional single RMSE over all values; it is a **sum of per-haplotype absolute errors** (since `sqrt((a-b)^2) = |a-b|`) averaged across events.

------

## Output

A single CSV summary is written to `self.output_dir`:

```
os.path.join(self.output_dir, outfile)
```

Default:

```
{self.output_dir}/mirror_subclone_result.csv
```

------

## Output Table Schema

One row per tool:

| Column | Meaning                                                      |
| ------ | ------------------------------------------------------------ |
| `Tool` | Tool name                                                    |
| `RMSE` | Mean aggregated haplotype absolute error across Clone1/Clone2 (as defined above) |
| `ACC`  | Exact-match accuracy requiring all four haplotype values to match |

------

## Return Value

Returns a `pd.DataFrame` with columns:

- `Tool`
- `RMSE`
- `ACC`

------

## Example

```python
from hcbench.gtbench.gtbench import GTBench

bench = GTBench(output_dir="out/gt_output")

df = bench.mirrorsubclone(
    tool_cna_files=[
        "/path/to/chisel/haplotype_combined.csv",
        "/path/to/signals/haplotype_combined.csv",
    ],
    tool_names=["CHISEL", "SIGNALS"],
    changes_file="/path/to/gt/mirrored_subclone_events.csv",
    profile_bin_size=100000,
    outfile="mirror_subclone_result.csv",
)

print(df)
```

