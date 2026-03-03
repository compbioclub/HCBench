## Function

```python
segmentation(
    self,
    gt_cna_file: str,
    tool_cna_files: List[str],
    tool_names: List[str],
    profile_bin_size=100000,
    threshold=0.95,
    outprefix: str = "segmentation_",
    level="cell",
    changes_file: Optional[str] = ""
)
```

This function evaluates **segmentation quality** of CNA calls by comparing tool predictions against **GT-derived CNA segments**.

It first derives CNA segments from the ground-truth CNA matrix using `annotate_segments(...)` (with optional CNV type detection). Then, for each tool, it converts the tool CNA matrix into a long format (region × cell), re-bins regions to a uniform resolution, and computes segmentation metrics using `get_segment_metric(...)`.

A single metrics CSV is written to disk.

------

## Parameters

| Name               | Type            | Description                                                  |
| ------------------ | --------------- | ------------------------------------------------------------ |
| `gt_cna_file`      | `str`           | Path to GT CNA matrix CSV (regions × cells). The first column is treated as index when reading. |
| `tool_cna_files`   | `List[str]`     | List of tool CNA profile CSVs. Must align with `tool_names`. |
| `tool_names`       | `List[str]`     | Tool names used to label metric outputs.                     |
| `profile_bin_size` | `int`           | Bin size used to split/standardize tool regions before metric computation. Default: `100000` (100kb). |
| `threshold`        | `float`         | Threshold passed into `annotate_segments(...)` (used for segmentation/CNV detection sensitivity). Default: `0.95`. |
| `outprefix`        | `str`           | Prefix used for the output metrics filename. Default: `"segmentation_"`. |
| `level`            | `str`           | Placeholder for evaluation level (`"cell"` or `"clone"`). Currently only `"cell"` path is active. |
| `changes_file`     | `Optional[str]` | Reserved for clone-level evaluation (commented out in current implementation). Default: `""`. |

------

## Input File Format

### `gt_cna_file` (GT CNA matrix)

Loaded as:

```python
gt_profile = pd.read_csv(gt_cna_file, index_col=0)
```

So the CSV is expected to be:

- rows: genomic regions (in the index)
- columns: cell IDs
- values: CNA states

Example:

```csv
region,cell_001,cell_002
chr1:1-100000,1|1,1|1
chr1:100001-200000,2|1,2|1
```

The GT matrix is transposed before segmentation annotation:

```python
gt_annotated_df = annotate_segments(gt_profile.T, detect_cnv_type=True, threshold=threshold)
```

`annotate_segments(...)` is expected to return a dataframe that includes (at minimum):

- `Chrom`
- `Start`
- `End`

A region string is then constructed:

```python
gt_annotated_df["region"] = (
    gt_annotated_df["Chrom"] + ":" +
    gt_annotated_df["Start"].astype(str) + "-" +
    gt_annotated_df["End"].astype(str)
)
```

### Tool CNA files (`tool_cna_files`)

Loaded via:

```python
tool_cna_df = read_and_drop_empty(f)
```

Expected to include a `region` column plus one column per cell (wide CNA matrix).

Example:

```csv
region,cell_001,cell_002
chr1:1-100000,1|1,1|1
chr1:100001-200000,2|1,2|1
```

------

## Processing Steps

For each tool:

1. Convert wide CNA matrix to long format (one row per region × cell):

```python
tool_cna_df_long = (
    pd.DataFrame(tool_cna_df)
      .melt(id_vars="region", var_name="cell", value_name="value")
      .dropna(subset=["value"])
)
```

1. Re-bin/split regions to `profile_bin_size` resolution:

```python
tool_cna_df_long = split_all_regions(tool_cna_df_long.set_index("region"), profile_bin_size)
tool_cna_df_long = tool_cna_df_long.reset_index().rename(columns={"index": "region"})
```

1. Compute segmentation metrics comparing tool calls to GT segments:

```python
result2 = get_segment_metric(gt_annotated_df, tool_cna_df_long)
result2["Tool"] = name
```

1. Collect per-tool metric outputs and concatenate:

```python
metric_df = pd.concat(metric_results, ignore_index=True)
```

------

## Output

A single CSV is written to `self.output_dir`:

```
os.path.join(self.output_dir, f"{outprefix}metrics.csv")
```

Default:

```
{self.output_dir}/segmentation_metrics.csv
```

The output rows depend on what `get_segment_metric(...)` returns (e.g., overlap scores, breakpoint precision/recall, CNV-type metrics, etc.), with an added column:

- `Tool`

------

## Return Value

This function does **not** explicitly return a value (runtime return is `None`).

------

## Example

```python
from hcbench.gtbench.gtbench import GTBench

bench = GTBench(output_dir="out/gt_output")

bench.segmentation(
    gt_cna_file="/path/to/gt/haplotype_combined.csv",
    tool_cna_files=[
        "/path/to/chisel/haplotype_combined.csv",
        "/path/to/signals/haplotype_combined.csv",
    ],
    tool_names=["CHISEL", "SIGNALS"],
    profile_bin_size=100000,
    threshold=0.8,
    outprefix="segmentation_",
)
```

