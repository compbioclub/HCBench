## Function

```python
cndetect(
    self,
    tool_cna_files: List[str],
    cna_profile_file: str,
    tool_names: List[str],
    haplotype: str = "combined",
    profile_bin_size=100000,
    outfile: str = "bin_level_results.csv",
    index_col: Optional[int] = 0,
) -> pd.DataFrame
```

This function evaluates **bin-level CNA detection performance** of multiple tools against a reference CNA profile file (ground truth).

For each tool CNA matrix, it:

1. Loads the tool CNA file and the reference CNA file.
2. Splits tool regions into a uniform binning resolution (`profile_bin_size`).
3. Aligns tool predictions and reference profiles on common regions and cells.
4. Computes evaluation metrics (**RMSE**, **SCC**, **ACC**) under a specified haplotype mode.
5. Saves a summary CSV across all tools.

------

## Parameters

| Name               | Type            | Description                                                  |
| ------------------ | --------------- | ------------------------------------------------------------ |
| `tool_cna_files`   | `List[str]`     | List of tool CNA profile CSV files. Must align with `tool_names` in order. |
| `cna_profile_file` | `str`           | Path to the reference (truth) CNA profile CSV file.          |
| `tool_names`       | `List[str]`     | Tool names used in logs and the result table.                |
| `haplotype`        | `str`           | Haplotype evaluation mode passed to `evaluate_haplotype_predictions`. Default: `"combined"`. |
| `profile_bin_size` | `int`           | Bin size used to split regions in tool predictions before alignment. Default: `100000` (100kb). |
| `outfile`          | `str`           | Output filename for the summary table. Default: `"bin_level_results.csv"`. |
| `index_col`        | `Optional[int]` | (Currently unused in implementation.) Intended CSV index column specification. Default: `0`. |

------

## Input File Format

### `cna_profile_file` (truth) and each file in `tool_cna_files`

Expected to be CNA matrices with:

- one column named `region`
- remaining columns as cell IDs
- values as CNA states (e.g., `"1|1"`, `"2|1"`), consistent with your pipeline conventions

Example:

```csv
region,cell_001,cell_002,cell_003
chr1:1-100000,1|1,1|1,2|1
chr1:100001-200000,1|1,1|1,2|1
```

Implementation details:

- Both truth and predictions are loaded with `read_and_drop_empty(...)` (drops empty columns/cells).
- The truth dataframe is indexed by `"region"` (`truth_df.set_index("region", inplace=True)`).
- Tool predictions are re-binned using:

```python
pred = split_all_regions(pred.set_index("region"), profile_bin_size)
```

Then truth and pred are aligned with:

```python
truth, pred = align(truth_df, pred)
```

------

## Evaluation Metrics

For each tool, metrics are computed via:

```python
rmse, scc, acc = evaluate_haplotype_predictions(pred, truth, haplotype)
```

The output table contains:

| Column | Meaning                                                      |
| ------ | ------------------------------------------------------------ |
| `Tool` | Tool name                                                    |
| `RMSE` | Root Mean Squared Error between prediction and truth         |
| `SCC`  | Similarity/Correlation metric returned by your evaluator (often “Spearman Correlation Coefficient” in CNA contexts, but exact definition depends on your implementation) |
| `ACC`  | Accuracy metric returned by your evaluator                   |

Notes:

- `SCC` may be `None` for certain inputs depending on your `evaluate_haplotype_predictions` implementation (e.g., constant vectors / insufficient comparable bins).

------

## Output

A single CSV is written to `self.output_dir`:

```
os.path.join(self.output_dir, outfile)
```

Default path:

```
{self.output_dir}/bin_level_results.csv
```

Example output:

```csv
Tool,RMSE,SCC,ACC
CHISEL,0.12,0.85,0.93
SIGNALS,0.20,0.71,0.88
```

------

## Return Value

Returns a `pd.DataFrame` with one row per tool and the columns:

- `Tool`
- `RMSE`
- `SCC`
- `ACC`

------

## Example

```python
from hcbench.gtbench.gtbench import GTBench

bench = GTBench(output_dir="out/gt_output")

df = bench.cndetect(
    tool_cna_files=[
        "/path/to/chisel/haplotype_combined.csv",
        "/path/to/signals/haplotype_combined.csv",
    ],
    cna_profile_file="/path/to/gt/haplotype_combined.csv",
    tool_names=["CHISEL", "SIGNALS"],
    haplotype="combined",
    profile_bin_size=100000,
    outfile="bin_level_results.csv",
)
print(df)
```

